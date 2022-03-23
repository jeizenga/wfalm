//
//  wf_test.cpp
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

#include <iostream>
#include <cmath>
#include <cassert>

#include "wfa_lm.hpp"
#include "wfa_lm_st.hpp"


template<int N>
int score_cigar(const std::string& seq1, const std::string& seq2,
                const std::vector<wfalm::CIGAROp>& cigar,
                int mismatch, std::array<uint32_t, N> gap_open, std::array<uint32_t, std::max(1, N)> gap_extend) {
    int s = 0;
    int i = 0;
    int j = 0;
    for (auto c : cigar) {
        switch (c.op) {
            case 'I':
            {
                int min_penalty = 100000000;
                for (size_t k = 0; k < gap_extend.size(); ++k) {
                    int penalty;
                    if (k < gap_open.size()) {
                        penalty = gap_open[k] + c.len * gap_extend[k];
                    }
                    else {
                        penalty = c.len * gap_extend[k];
                    }
                    min_penalty = std::min(min_penalty, penalty);
                }
                i += c.len;
                s += min_penalty;
                break;
            }
            case 'D':
            {
                int min_penalty = 100000000;
                for (size_t k = 0; k < gap_extend.size(); ++k) {
                    int penalty;
                    if (k < gap_open.size()) {
                        penalty = gap_open[k] + c.len * gap_extend[k];
                    }
                    else {
                        penalty = c.len * gap_extend[k];
                    }
                    min_penalty = std::min(min_penalty, penalty);
                }
                j += c.len;
                s += min_penalty;
                break;
            }
            case 'M':
            {
                for (int k = 0; k < c.len; ++k, ++i, ++j) {
                    if (seq1[i] != seq2[j]) {
                        s += mismatch;
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
    
    // anchor            |       |
    //  0    2   3    4  6        7  8     11     12
    //           *    *           *  *            *
    //  TACGGTCACCGCGACGA-CGACGGCAATTACTCCAAGTTGTCT
    //  TACGG-CACGGCGATGAACGACGGCACTTTCTCCA--TTGTCC
    //                   ^2 placements for this
    // 0     1   1    2  3        5  5     67     8
    //       1   9    9  4        2  8     90     2
    //      1   1    2  3        5  5     66     88
    //      0   7    7  3        0  6     89     02
    std::string seq1 = "TACGGTCACCGCGACGACGACGGCAATTACTCCAAGTTGTCT";
    std::string seq2 = "TACGGCACGGCGATGAACGACGGCACTTTCTCCATTGTCC";

//    uint32_t mismatch = 1;
//    std::array<uint32_t, 0> gap_open;
//    std::array<uint32_t, 1> gap_extend{1};
//    uint32_t max_mem = 400;
    
    uint32_t mismatch = 1;
    std::array<uint32_t, 1> gap_open{1};
    std::array<uint32_t, 1> gap_extend{1};
    uint32_t max_mem = 400;

//    // TACGGC-CGGCGATG--------CACTTTCTCCATTG-CC
//    // TACGGCACGGCGATGAACGACGGCACTTTCTCCATTGTCC
//    std::string seq1 = "TACGGCCGGCGATGCACTTTCTCCATTGCC";
//    std::string seq2 = "TACGGCACGGCGATGAACGACGGCACTTTCTCCATTGTCC";
//
//    uint32_t mismatch = 1;
//    std::array<uint32_t, 2> gap_open{1, 3};
//    std::array<uint32_t, 2> gap_extend{2, 1};
//    uint32_t max_mem = 400;

    wfalm::WFAligner<1> aligner(mismatch, gap_extend, gap_open);
    wfalm::WFAlignerST<1> aligner_st(mismatch, gap_extend, gap_open);
    
    int prune = -1;//ceil(-2.0 * log(1e-12) / log(4.0));
    aligner.lagging_diagonal_prune = prune;
    
    std::cout << "GLOBAL ALIGNMENT" << std::endl;

    std::cout << "sequence 1: " << std::endl;
    std::cout << seq1 << std::endl;
    std::cout << "sequence 2: " << std::endl;
    std::cout << seq2 << std::endl;
    std::cout << "params:" << std::endl;
    std::cout << "\tx = " << mismatch << std::endl;
    std::cout << "\to =";
    for (auto o : gap_open) {
        std::cout << ' ' << o;
    }
    std::cout << std::endl;
    std::cout << "\te =";
    for (auto e : gap_extend) {
        std::cout << ' ' << e;
    }
    std::cout << std::endl;
    std::cout << "\tp = " << prune << std::endl;
    std::cout << "\tM = " << max_mem << std::endl;
    for (int mem : {0, 1, 2, 3}) {
        for (bool use_st : {false, true}) {
            bool low_mem = (mem == 1);
            bool recursive = (mem == 2);
            bool adaptive = (mem == 3);
            std::pair<std::vector<wfalm::CIGAROp>, int32_t> res;
            if (recursive && use_st) {
                std::cout << "recursive algorithm with suffix tree match:" << std::endl;
                res = aligner_st.wavefront_align_recursive(seq1.c_str(), seq1.size(),
                                                           seq2.c_str(), seq2.size());
            }
            else if (recursive) {
                std::cout << "recursive algorithm:" << std::endl;
                res = aligner.wavefront_align_recursive(seq1.c_str(), seq1.size(),
                                                        seq2.c_str(), seq2.size());
            }
            else if (adaptive && use_st) {
                std::cout << "adaptive algorithm with suffix tree match:" << std::endl;
                res = aligner_st.wavefront_align_adaptive(seq1.c_str(), seq1.size(),
                                                          seq2.c_str(), seq2.size(),
                                                          max_mem);
            }
            else if (adaptive) {
                std::cout << "adaptive algorithm:" << std::endl;
                res = aligner.wavefront_align_adaptive(seq1.c_str(), seq1.size(),
                                                       seq2.c_str(), seq2.size(),
                                                       max_mem);
            }
            else if (low_mem && use_st) {
                std::cout << "low memory algorithm with suffix tree match:" << std::endl;
                res = aligner_st.wavefront_align_low_mem(seq1.c_str(), seq1.size(),
                                                         seq2.c_str(), seq2.size());
            }
            else if (low_mem) {
                std::cout << "low memory algorithm:" << std::endl;
                res = aligner.wavefront_align_low_mem(seq1.c_str(), seq1.size(),
                                                      seq2.c_str(), seq2.size());
            }
            else if (use_st) {
                std::cout << "standard algorithm with suffix tree match:" << std::endl;
                res = aligner_st.wavefront_align(seq1.c_str(), seq1.size(),
                                                 seq2.c_str(), seq2.size());
            }
            else {
                std::cout << "standard algorithm:" << std::endl;
                res = aligner.wavefront_align(seq1.c_str(), seq1.size(),
                                              seq2.c_str(), seq2.size());

            }
            auto cigar = res.first;
            int s = res.second;

            std::cout << "score: " << s << std::endl;
            std::cout << "score density: " << double(s) / double(seq1.size() + seq2.size()) << std::endl;
            std::cout << "CIGAR string:" << std::endl;
            for (auto cigar_op : cigar) {
                std::cout << cigar_op.len << cigar_op.op;
            }
            std::cout << std::endl;
            
            assert(s == score_cigar<1>(seq1, seq2, cigar, mismatch, gap_open, gap_extend));
        }
    }
    
    std::string pref1 = "CACCATAAGTCAGG";
    std::string pref2 = "ACCCGGCTCACGCC";
    std::string suff1 = "GAGAACTCTAATGG";
    std::string suff2 = "CTTTTAAATGGTAG";

    auto local_seq1 = pref1 + seq1 + suff1;
    auto local_seq2 = pref2 + seq2 + suff2;

    uint32_t match_sw = 1;
    uint32_t mismatch_sw = 1;
    uint32_t gap_open_sw = 2;
    uint32_t gap_extend_sw = 1;
    int max_mem_sw = 200;

    wfalm::WFAligner<1> aligner_sw(match_sw, mismatch_sw, {gap_extend_sw}, {gap_open_sw});
    wfalm::WFAlignerST<1> aligner_sw_st(match_sw, mismatch_sw, {gap_extend_sw}, {gap_open_sw});
    
//    uint32_t match_sw = 1;
//    uint32_t mismatch_sw = 1;
//    uint32_t gap_extend_sw = 1;
//    uint32_t gap_open_sw = 0;
//    int max_mem_sw = 200;
//
//    wfalm::WFAligner<0> aligner_sw(match_sw, mismatch_sw, {gap_extend_sw}, {});
//    wfalm::WFAlignerST<0> aligner_sw_st(match_sw, mismatch_sw, {gap_extend_sw}, {});

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
    std::cout << "\tM = " << max_mem_sw << std::endl;
    for (int mem : {0, 1, 2, 3}) {
        for (bool use_st : {false, true}) {
            bool low_mem = (mem == 1);
            bool recursive = (mem == 2);
            bool adaptive = (mem == 3);
            std::tuple<std::vector<wfalm::CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>> res;
            if (low_mem && !use_st) {
                std::cout << "low memory algorithm:" << std::endl;
                res = aligner_sw.wavefront_align_local_low_mem(local_seq1.c_str(), local_seq1.size(),
                                                               local_seq2.c_str(), local_seq2.size(),
                                                               anchor_begin_1, anchor_end_1,
                                                               anchor_begin_2, anchor_end_2,
                                                               false);
            }
            else if (low_mem && use_st) {
                std::cout << "low memory algorithm with suffix tree match:" << std::endl;
                res = aligner_sw_st.wavefront_align_local_low_mem(local_seq1.c_str(), local_seq1.size(),
                                                                  local_seq2.c_str(), local_seq2.size(),
                                                                  anchor_begin_1, anchor_end_1,
                                                                  anchor_begin_2, anchor_end_2,
                                                                  false);
            }
            else if (adaptive && use_st) {
                std::cout << "adaptive algorithm with suffix tree match:" << std::endl;
                res = aligner_sw_st.wavefront_align_local_adaptive(local_seq1.c_str(), local_seq1.size(),
                                                                   local_seq2.c_str(), local_seq2.size(),
                                                                   max_mem_sw,
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   false);
            }
            else if (adaptive && !use_st) {
                std::cout << "adaptive algorithm:" << std::endl;
                res = aligner_sw.wavefront_align_local_adaptive(local_seq1.c_str(), local_seq1.size(),
                                                                local_seq2.c_str(), local_seq2.size(),
                                                                max_mem_sw,
                                                                anchor_begin_1, anchor_end_1,
                                                                anchor_begin_2, anchor_end_2,
                                                                false);
            }
            else if (recursive && use_st) {
                std::cout << "recursive algorithm with suffix tree match:" << std::endl;
                res = aligner_sw_st.wavefront_align_local_recursive(local_seq1.c_str(), local_seq1.size(),
                                                                    local_seq2.c_str(), local_seq2.size(),
                                                                    anchor_begin_1, anchor_end_1,
                                                                    anchor_begin_2, anchor_end_2,
                                                                    false);
            }
            else if (recursive && !use_st) {
                std::cout << "recursive algorithm:" << std::endl;
                res = aligner_sw.wavefront_align_local_recursive(local_seq1.c_str(), local_seq1.size(),
                                                                 local_seq2.c_str(), local_seq2.size(),
                                                                 anchor_begin_1, anchor_end_1,
                                                                 anchor_begin_2, anchor_end_2,
                                                                 false);
            }
            else if (use_st) {
                std::cout << "standard algorithm with suffix tree match:" << std::endl;
                res = aligner_sw_st.wavefront_align_local(local_seq1.c_str(), local_seq1.size(),
                                                          local_seq2.c_str(), local_seq2.size(),
                                                          anchor_begin_1, anchor_end_1,
                                                          anchor_begin_2, anchor_end_2,
                                                          false);
            }
            else {
                std::cout << "standard algorithm:" << std::endl;
                res = aligner_sw.wavefront_align_local(local_seq1.c_str(), local_seq1.size(),
                                                       local_seq2.c_str(), local_seq2.size(),
                                                       anchor_begin_1, anchor_end_1,
                                                       anchor_begin_2, anchor_end_2,
                                                       false);
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
