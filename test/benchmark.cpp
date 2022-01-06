//
//  benchmark.cpp
//  
// MIT License
//
// Copyright (c) 2022 Jordan Eizenga
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <chrono>
#include <cassert>
#include <sys/resource.h>

#include "wfa_lm.hpp"
#include "wfa_lm_st.hpp"

using namespace std;
using namespace std::chrono;
using namespace wfalm;

int default_match = 0;
int default_mismatch = 4;
int default_gap_open = 5;
int default_gap_extend = 1;

void print_usage() {
    cerr << "usage: benchmark [options] seq1.fa seq2.fa" << endl;
    cerr << endl;
    cerr << "options:" << endl;
    cerr << " -s, --standard-mem    O(s^2) memory algorithm [default]" << endl;
    cerr << " -m, --low-mem         O(s^3/2) memory algrorithm" << endl;
    cerr << " -g, --global          global alignment [default]" << endl;
    cerr << " -l, --local SEED      local alignment, requires seed of format i:j,k:l" << endl;
    cerr << " -c, --compare-match   compute matches by direct comparison, O(sN) time [default]" << endl;
    cerr << " -t, --suf-tree-match  compute matches with suffix tree, O(s^2 + N) time" << endl;
    cerr << " -x, --mismatch INT    mismatch score [default " << default_mismatch << "] " << endl;
    cerr << " -o, --gap-open INT    gap open score [default " << default_mismatch << "] " << endl;
    cerr << " -e, --gap-extend INT  gap extend score [default " << default_mismatch << "] " << endl;
    cerr << " -a, --match INT       match score (triggers SWG-style params) [default " << default_mismatch << "] " << endl;
    cerr << " -h, --help            print this message and exit" << endl;
}


void parse_seed_string(const string& str, size_t& anchor11, size_t& anchor12,
                       size_t& anchor21, size_t& anchor22) {
    
    auto sep = str.find(",");
    if (sep == string::npos) {
        cerr << "ERROR: malformed seed string" << endl;
        exit(1);
    }
    
    auto seed1 = str.substr(0, sep);
    auto seed2 = str.substr(sep + 1, string::npos);
    
    auto colon1 = seed1.find(":");
    auto colon2 = seed2.find(":");
    
    if (colon1 == string::npos || colon2 == string::npos) {
        cerr << "ERROR: malformed seed string" << endl;
        exit(1);
    }
    
    anchor11 = atoi(seed1.substr(0, colon1).c_str());
    anchor12 = atoi(seed1.substr(colon1 + 1, string::npos).c_str());
    anchor21 = atoi(seed2.substr(0, colon2).c_str());
    anchor22 = atoi(seed2.substr(colon2 + 1, string::npos).c_str());
}

string cigar_to_string(const vector<CIGAROp>& cigar) {
    stringstream strm;
    for (auto op : cigar) {
        strm << op.len << op.op;
    }
    return strm.str();
}

string parse_fasta(const string& fasta_name, string& seq_name) {
    
    string seq;
    
    string line;
    
    ifstream strm(fasta_name);
    bool seen_header = false;
    while (strm) {
        line.clear();
        getline(strm, line);
        if (!line.empty() && line.front() == '>') {
            if (!seen_header) {
                seq_name = line.substr(1, line.find(' ') - 1);
                cerr << "loading sequence " << seq_name << " from " << fasta_name << endl;
                seen_header = true;
            }
            else {
                cerr << "WARNING: FASTA file contains multiple sequences. Only the first one will be aligned." << endl;
                break;
            }
        }
        else {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq.append(line);
        }
    }
    
    return seq;
}


string shortened_seq(const string& seq, size_t len) {
    string shortened;
    if (seq.size() <= len) {
        shortened = seq;
    }
    else {
        size_t l1 = len / 2;
        size_t l2 = len - l1;
        shortened = seq.substr(0, l1) + string("...") + seq.substr(seq.size() - l2, l2);
    }
    return shortened;
}

int64_t max_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return 1024 * usage.ru_maxrss;
}

int main(int argc, char **argv) {
    
    int match = default_match;
    int mismatch = default_mismatch;
    int gap_open = default_gap_open;
    int gap_extend = default_gap_extend;
    size_t anchor11 = 0, anchor12 = 0, anchor21 = 0, anchor22 = 0;
    bool local = false;
    bool compare_match = true;
    bool low_memory = false;
    
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"standard-mem", no_argument, 0, 's'},
            {"low-mem", no_argument, 0, 'm'},
            {"global", no_argument, 0, 'g'},
            {"local", required_argument, 0, 'l'},
            {"compare-match", no_argument, 0, 'c'},
            {"suf-tree-match", no_argument, 0, 't'},
            {"mismatch", required_argument, 0, 'x'},
            {"gap-open", required_argument, 0, 'o'},
            {"gap-extend", required_argument, 0, 'e'},
            {"match", required_argument, 0, 'a'},
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int option_index = 0;
        c = getopt_long (argc, argv, "sml:ctx:o:e:a:h",
                         long_options, &option_index);
        if (c == -1){
            break;
        }
        
        switch(c){
            case 's':
                low_memory = false;
                break;
            case 'm':
                low_memory = true;
                break;
            case 'g':
                local = false;
                break;
            case 'l':
                local = true;
                parse_seed_string(optarg, anchor11, anchor12, anchor21, anchor22);
                break;
            case 'c':
                compare_match = true;
                break;
            case 't':
                compare_match = false;
                break;
            case 'x':
                mismatch = atoi(optarg);
                break;
            case 'o':
                gap_open = atoi(optarg);
                break;
            case 'e':
                gap_extend = atoi(optarg);
                break;
            case 'a':
                match = atoi(optarg);
                break;
            case 'h':
            case '?':
                print_usage();
                return 0;
            default:
                print_usage();
                return 1;
        }
    }
    
    if (optind + 2 != argc) {
        // no positional arguments
        cerr << "ERROR: expected 2 position arguments, but got " << (argc - optind) << endl;
        print_usage();
        return 1;
    }
    if (match < 0 || mismatch < 0 || gap_open < 0 || gap_extend < 0) {
        cerr << "ERROR: all scoring params must be non-negative" << endl;
        return 1;
    }
    if (match == 0 && local) {
        cerr << "ERROR: must use Smith-Waterman-Gotoh style parameters in local alignment" << endl;
        return 1;
    }
    
    string fasta1(argv[optind++]);
    string fasta2(argv[optind++]);
    
    cerr << "options: " << (low_memory ? "low" : "standard") << " memory, " << (compare_match ? "direct" : "suffix tree") << " match function, " << (local ? "local" : "global") << " alignment" << endl;
    if (match != 0) {
        cerr << "a = " << match << ", ";
    }
    cerr << "x = " << mismatch << ", o = " << gap_open << ", e = " << gap_extend << endl;
    cerr << "reading FASTA files " << fasta1 << " and " << fasta2 << endl;
    
    steady_clock::time_point begin_parse = steady_clock::now();
    
    string seqname1, seqname2;
    string seq1 = parse_fasta(fasta1, seqname1);
    string seq2 = parse_fasta(fasta2, seqname2);
    
    steady_clock::time_point end_parse = steady_clock::now();
    
    auto parse_time = duration_cast<microseconds>(end_parse - begin_parse);
    cerr << "parsing completed in " << parse_time.count() << " us" << endl;
    //cerr << "memory usage " << max_memory_usage() << endl;
    
    cerr << "aligning sequences '" << seqname1 << "' and '" << seqname2 << "':" << endl;
    cerr << shortened_seq(seq1, 50) << endl;
    cerr << shortened_seq(seq2, 50) << endl;
    
    
    WFScores wf_scores(mismatch, gap_open, gap_extend);
    SWGScores swg_scores(match, mismatch, gap_open, gap_extend);
    
    vector<CIGAROp> cigar;
    int32_t score;
    pair<size_t, size_t> range1;
    pair<size_t, size_t> range2;
    
    steady_clock::time_point begin_align = steady_clock::now();
    
    if (!low_memory && !local && compare_match) {
        if (match == 0) {
            tie(cigar, score) = wavefront_align(seq1, seq2, wf_scores);
        }
        else {
            tie(cigar, score) = wavefront_align(seq1, seq2, swg_scores);
        }
    }
    else if (low_memory && !local && compare_match) {
        if (match == 0) {
            tie(cigar, score) = wavefront_align_low_mem(seq1, seq2, wf_scores);
        }
        else {
            tie(cigar, score) = wavefront_align_low_mem(seq1, seq2, swg_scores);
        }
    }
    else if (!low_memory && !local && !compare_match) {
        if (match == 0) {
            tie(cigar, score) = wavefront_align_st(seq1, seq2, wf_scores);
        }
        else {
            tie(cigar, score) = wavefront_align_st(seq1, seq2, swg_scores);
        }
    }
    else if (low_memory && !local && !compare_match) {
        if (match == 0) {
            tie(cigar, score) = wavefront_align_low_mem_st(seq1, seq2, wf_scores);
        }
        else {
            tie(cigar, score) = wavefront_align_low_mem_st(seq1, seq2, swg_scores);
        }
    }
    else if (!low_memory && local && compare_match) {
        tie(cigar, score, range1, range2) = wavefront_align_local(seq1, seq2, anchor11, anchor12, anchor21, anchor22, swg_scores);
    }
    else if (low_memory && local && compare_match) {
        tie(cigar, score, range1, range2) = wavefront_align_local_low_mem(seq1, seq2, anchor11, anchor12, anchor21, anchor22, swg_scores);
    }
    else if (!low_memory && local && !compare_match) {
        tie(cigar, score, range1, range2) = wavefront_align_local_st(seq1, seq2, anchor11, anchor12, anchor21, anchor22, swg_scores);
    }
    else if (low_memory && local && !compare_match) {
        tie(cigar, score, range1, range2) = wavefront_align_local_low_mem_st(seq1, seq2, anchor11, anchor12, anchor21, anchor22, swg_scores);
    }
    
    steady_clock::time_point end_align = steady_clock::now();
    auto align_time = duration_cast<microseconds>(end_align - begin_align);
    
    cerr << "alignment:" << endl;
    cerr << shortened_seq(cigar_to_string(cigar), 400) << endl;
    if (local) {
        cerr << "aligned intervals: [" << range1.first << ":" << range1.second << "] [" << range2.first << ":" << range2.second << "]" << endl;
    }
    cerr << "alignment completed in " << align_time.count() << " us" << endl;
    //cerr << "memory usage " << max_memory_usage() << endl;
}
