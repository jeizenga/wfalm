//
//  wf_test.cpp
//  

#include <iostream>
#include <cmath>
#include <cassert>

#include "wfa_lm.hpp"
#include "wfa_lm_st.hpp"

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
    
    std::cout << "GLOBAL ALIGNMENT" << std::endl;

    std::cout << "sequence 1: " << std::endl;
    std::cout << seq1 << std::endl;
    std::cout << "sequence 2: " << std::endl;
    std::cout << seq2 << std::endl;
    std::cout << "params:" << std::endl;
    std::cout << "\tx = " << mismatch << std::endl;
    std::cout << "\to = " << gap_open << std::endl;
    std::cout << "\te = " << gap_extend << std::endl;
    std::cout << "\tp = " << prune << std::endl;
    for (bool low_mem : {false, true}) {
        for (bool use_st : {false, true}) {
            std::pair<std::vector<wfalm::CIGAROp>, int32_t> res;
            if (low_mem && use_st) {
                std::cout << "low memory algorithm with suffix tree match:" << std::endl;
                res = wfalm::wavefront_align_low_mem_st(seq1, seq2, scores, prune);
            }
            else if (low_mem) {
                std::cout << "low memory algorithm:" << std::endl;
                res = wfalm::wavefront_align_low_mem(seq1, seq2, scores, prune);
            }
            else if (use_st) {
                std::cout << "standard algorithm with suffix tree match:" << std::endl;
                res = wfalm::wavefront_align_st(seq1, seq2, scores, prune);
            }
            else {
                std::cout << "standard algorithm:" << std::endl;
                res = wfalm::wavefront_align(seq1, seq2, scores, prune);
                
            }
            auto cigar = res.first;
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
    }
    
    std::string pref1 = "CACCATAAGTCAGG";
    std::string pref2 = "ACCCGGCTCACGCC";
    std::string suff1 = "GAGAACTCTAATGG";
    std::string suff2 = "CTTTTAAATGGTAG";
    
    auto local_seq1 = pref1 + seq1 + suff1;
    auto local_seq2 = pref2 + seq2 + suff2;
    
    int match_sw = 1;
    int mismatch_sw = 1;
    int gap_open_sw = 2;
    int gap_extend_sw = 1;
    
    wfalm::SWGScores sw_scores(match_sw, mismatch_sw, gap_open_sw, gap_extend_sw);
    
    int anchor_begin_1 = pref1.size() + 16;
    int anchor_end_1 = pref1.size() + 25;
    int anchor_begin_2 = pref2.size() + 16;
    int anchor_end_2 = pref2.size() + 25;
    
    std::cout << std::endl << "LOCAL ALIGNMENT" << std::endl;
    
    std::cout << "sequence 1: " << std::endl;
    std::cout << local_seq1 << std::endl;
    std::cout << "sequence 2: " << std::endl;
    std::cout << local_seq2 << std::endl;
    std::cout << "params:" << std::endl;
    std::cout << "\tm = " << match_sw << std::endl;
    std::cout << "\tx = " << mismatch_sw << std::endl;
    std::cout << "\to = " << gap_open_sw << std::endl;
    std::cout << "\te = " << gap_extend_sw << std::endl;
    for (bool low_mem : {false, true}) {
        for (bool use_st : {false, true}) {
            std::tuple<std::vector<wfalm::CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>> res;
            if (low_mem && use_st) {
                std::cout << "low memory algorithm with suffix tree match:" << std::endl;
                res = wfalm::wavefront_align_local_low_mem_st(local_seq1, local_seq2,
                                                              anchor_begin_1, anchor_end_1,
                                                              anchor_begin_2, anchor_end_2,
                                                              sw_scores, false);
            }
            else if (low_mem) {
                std::cout << "low memory algorithm:" << std::endl;
                res = wfalm::wavefront_align_local_low_mem(local_seq1, local_seq2,
                                                           anchor_begin_1, anchor_end_1,
                                                           anchor_begin_2, anchor_end_2,
                                                           sw_scores, false);
            }
            else if (use_st) {
                std::cout << "standard algorithm with suffix tree match:" << std::endl;
                res = wfalm::wavefront_align_local_st(local_seq1, local_seq2,
                                                      anchor_begin_1, anchor_end_1,
                                                      anchor_begin_2, anchor_end_2,
                                                      sw_scores, false);
            }
            else {
                std::cout << "standard algorithm:" << std::endl;
                res = wfalm::wavefront_align_local(local_seq1, local_seq2,
                                                   anchor_begin_1, anchor_end_1,
                                                   anchor_begin_2, anchor_end_2,
                                                   sw_scores, false);
            }
            auto cigar = std::get<0>(res);
            
            std::cout << "score: " << std::get<1>(res) << std::endl;
            std::cout << "aligned intervals: s1[" << std::get<2>(res).first << ":" << std::get<2>(res).second << "], s2[" << std::get<3>(res).first << ":" << std::get<3>(res).second << "]" << std::endl;
            std::cout << "CIGAR string:" << std::endl;
            for (auto cigar_op : cigar) {
                std::cout << cigar_op.len << cigar_op.op;
            }
            std::cout << std::endl;
        }
    }
}
