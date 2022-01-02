//
//  wf_test.cpp
//  

#include <iostream>
#include <cmath>
#include <cassert>

#include "wfa_lm.hpp"

int score_cigar(const std::string& seq1, const std::string& seq2,
                const std::vector<wfalm::CIGAROp>& cigar, const wfalm::WFScores& scores) {
    int s = 0;
    int i = 0;
    int j = 0;
    for (auto c : cigar) {
        switch (c.op) {
            case 'I':
                i += c.len;
                s += scores.gap_open + c.len * scores.gap_extend;
                break;
            case 'D':
                j += c.len;
                s += scores.gap_open + c.len * scores.gap_extend;
                break;
            case 'M':
            {
                for (int k = 0; k < c.len; ++k, ++i, ++j) {
                    if (seq1[i] != seq2[j]) {
                        s += scores.mismatch;
                    }
                }
            }
                
            default:
                break;
        }
    }
    assert(cigar.empty() || i == seq1.size());
    assert(cigar.empty() || j == seq2.size());
    return s;
}

int main(int argc, char** argv) {
    
    // 0    2   3   4  6         7  8     11     12
    //          *   *            *  *            *
    // TACGGTCACCGCGACG-ACGACGGCAATTACTCCAAGTTGTCT
    // TACGG-CACGGCGATGAACGACGGCACTTTCTCCA--TTGTCC
    //                 ^2 placements for this
    std::string seq1 = "TACGGTCACCGCGACGACGACGGCAATTACTCCAAGTTGTCT";
    std::string seq2 = "TACGGCACGGCGATGAACGACGGCACTTTCTCCATTGTCC";
    
    int mismatch = 1;
    int gap_open = 1;
    int gap_extend = 1;
    wfalm::WFScores scores(mismatch, gap_open, gap_extend);
    
    int prune = -1;//ceil(-2.0 * log(1e-12) / log(4.0));
    auto res = wfalm::wavefront_align_low_mem(seq1, seq2, scores, prune);
    std::vector<wfalm::CIGAROp> cigar = res.first;
    
    std::cout << "sequence 1: " << std::endl;
    std::cout << seq1 << std::endl;
    std::cout << "sequence 2: " << std::endl;
    std::cout << seq2 << std::endl;
    std::cout << "params:" << std::endl;
    std::cout << "\tx = " << mismatch << std::endl;
    std::cout << "\to = " << gap_open << std::endl;
    std::cout << "\te = " << gap_extend << std::endl;
    std::cout << "\tp = " << prune << std::endl;
    int s = score_cigar(seq1, seq2, cigar, scores);
    assert(s == res.second);
    std::cout << "score: " << s << std::endl;
    std::cout << "score density: " << double(s) / double(seq1.size() + seq2.size()) << std::endl;
    std::cout << "CIGAR string:" << std::endl;
    for (auto cigar_op : cigar) {
        std::cout << cigar_op.len << cigar_op.op;
    }
    std::cout << std::endl;
}
