//
//  wfa_lm.hpp
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


#ifndef wfa_lm_hpp_included
#define wfa_lm_hpp_included

#include <cstdint>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>
#include <deque>
#include <sstream>
#include <limits>

namespace wfalm {

/*
 * Score parameters for wavefront alignment. On opening a insertion or
 * deletion, *both* the gap open and gap extend penalties are applied.
 */
struct WFScores {
    WFScores(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    WFScores() = delete;
    
    uint32_t mismatch;
    uint32_t gap_open;
    uint32_t gap_extend;
};

/*
 * Score parameters for Smith-Waterman-Gotoh or Needleman-Wunsch alignment.
 * On opening a insertion or deletion, *both* the gap open and gap extend
 * penalties are applied.
 */
struct SWGScores {
    SWGScores(uint32_t match, uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    SWGScores() = delete;
    
    uint32_t match;
    uint32_t mismatch;
    uint32_t gap_open;
    uint32_t gap_extend;
};

/*
 * Operation in a CIGAR string (as defined in SAM format specification)
 */
struct CIGAROp {
    CIGAROp(char op, uint32_t len);
    CIGAROp() = default;
    
    char op;
    uint32_t len;
    
};

/// Align two sequences using the wavefront alignment algorithm, returns a CIGAR
/// string.
/// Optionally performs pruning (if prune_diff >= 0) of diagonals that are the
/// indicated difference behind the leading diagonal, measured in antidiagonals.
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align(const std::string& seq1, const std::string& seq2,
                const WFScores& scores, int32_t prune_diff = -1);

/// Same semantics as above, but with lower memory complexity
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem(const std::string& seq1, const std::string& seq2,
                        const WFScores& scores, int32_t prune_diff = -1);


/*********************************
 * Only internal functions below here
 *********************************/










WFScores::WFScores(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend)
    : mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend)
{
    
}

SWGScores::SWGScores(uint32_t match, uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend)
    : match(match), mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend)
{
    
}

CIGAROp::CIGAROp(char op, uint32_t len) : op(op), len(len) {
    
}

namespace internal {

enum WFMatrix_t {MAT_M, MAT_D, MAT_I};

struct WFEntry {
    WFEntry() :
        M(std::numeric_limits<int32_t>::min()),
        I(std::numeric_limits<int32_t>::min()),
        D(std::numeric_limits<int32_t>::min())
    {
        
    }
    int32_t M;
    int32_t I;
    int32_t D;
};

// TODO: replace deque with manual pointers? (will need to implement ctors)
struct Wavefront {
public:
    Wavefront() {}
    Wavefront(int32_t diag_begin, int32_t diag_end)
    : diag_begin(diag_begin), entries(diag_end - diag_begin) {
        
    }
    
    int32_t diag_begin;
    std::deque<WFEntry> entries;
};

inline
void wavefront_extend(const std::string& seq1, const std::string& seq2,
                      Wavefront& wf) {
    
    for (int64_t k = 0; k < wf.entries.size(); ++k) {
        int64_t diag = wf.diag_begin + k;
        int64_t anti_diag = wf.entries[k].M;
        int64_t i = (diag + anti_diag) / 2 + 1;
        int64_t j = (anti_diag - diag) / 2 + 1;
        while (i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]) {
            wf.entries[k].M += 2;
            ++i;
            ++j;
        }
    }
}

template<typename WFVector>
Wavefront wavefront_next(const std::string& seq1, const std::string& seq2,
                         const WFScores& scores, const WFVector& wfs) {
    
    // find the range of incoming wavefronts
    int32_t hi = std::numeric_limits<int32_t>::min();
    int32_t lo = std::numeric_limits<int32_t>::max();
    // mismatch stays in the same diagonals
    if (wfs.size() >= scores.mismatch) {
        const auto& wf_prev = wfs[wfs.size() - scores.mismatch];
        if (!wf_prev.entries.empty()) {
            lo = std::min<int32_t>(lo, wf_prev.diag_begin);
            hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.entries.size());
        }
    }
    // TODO: but sometimes these bounds are slightly too large because of the gap extend
    // insertion/deletion expand the diagonals by 1
    for (auto s_back : {scores.gap_extend, scores.gap_open + scores.gap_extend}) {
        if (wfs.size() >= s_back) {
            const auto& wf_prev = wfs[wfs.size() - s_back];
            if (!wf_prev.entries.empty()) {
                lo = std::min<int32_t>(lo, wf_prev.diag_begin - 1);
                hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.entries.size() + 1);
            }
        }
    }
    
    // some twiddly code so that we can use the same constructor and preserve RVO
    int32_t ctor_lo, ctor_hi;
    if (hi < lo) {
        ctor_lo = 0;
        ctor_hi = 0;
    }
    else {
        ctor_lo = lo;
        ctor_hi = hi;
    }
    
    Wavefront wf(ctor_lo, ctor_hi);
    
    if (lo <= hi) {
        
        // open gaps
        if (wfs.size() >= scores.gap_extend + scores.gap_open) {
            const auto& wf_prev = wfs[wfs.size() - scores.gap_extend - scores.gap_open];
            for (auto k = lo; k < hi; ++k) {
                if (k - 1 >= wf_prev.diag_begin && k - 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    int64_t a = wf_prev.entries[k - 1 - wf_prev.diag_begin].M + 1;
                    if ((a + k) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                        wf.entries[k - lo].I = a;
                    }
                }
                if (k + 1 >= wf_prev.diag_begin && k + 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    int64_t a = wf_prev.entries[k + 1 - wf_prev.diag_begin].M + 1;
                    if ((a + k) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                        wf.entries[k - lo].D = a;
                    }
                }
            }
        }
        
        // extend gaps
        if (wfs.size() >= scores.gap_extend) {
            const auto& wf_prev = wfs[wfs.size() - scores.gap_extend];
            for (auto k = lo; k < hi; ++k) {
                if (k - 1 >= wf_prev.diag_begin && k - 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    int64_t a = wf_prev.entries[k - 1 - wf_prev.diag_begin].I + 1;
                    if ((a + k) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                        wf.entries[k - lo].I = std::max<int32_t>(wf.entries[k - lo].I, a);
                    }
                }
                if (k + 1 >= wf_prev.diag_begin && k + 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    int64_t a = wf_prev.entries[k + 1 - wf_prev.diag_begin].D + 1;
                    if ((a + k) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                        wf.entries[k - lo].D = std::max<int32_t>(wf.entries[k - lo].D, a);
                    }
                }
            }
        }
        
        // mismatches
        if (wfs.size() >= scores.mismatch) {
            const auto& wf_prev = wfs[wfs.size() - scores.mismatch];
            for (auto k = lo; k < hi; ++k) {
                if (k >= wf_prev.diag_begin && k < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    int64_t a = wf_prev.entries[k - wf_prev.diag_begin].M + 2;
                    if ((a + k) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                        wf.entries[k - lo].M = a;
                    }
                }
            }
        }
        
        // calculate the opt
        for (size_t k = 0; k < wf.entries.size(); ++k) {
            wf.entries[k].M = std::max(wf.entries[k].M, std::max(wf.entries[k].I, wf.entries[k].D));
        }
    }
    
    return wf;
}

inline
void wavefront_prune(Wavefront& wf, int32_t sub_opt_diff) {
    
    if (!wf.entries.empty()) {
        int32_t max_anti_diag = std::numeric_limits<int32_t>::min();
        for (const auto& entry : wf.entries) {
            max_anti_diag = std::max(max_anti_diag, entry.M);
        }
        
        while (wf.entries.front().M < max_anti_diag - sub_opt_diff) {
            wf.entries.pop_front();
            ++wf.diag_begin;
        }
        
        while (wf.entries.back().M < max_anti_diag - sub_opt_diff) {
            wf.entries.pop_back();
        }
    }
}

inline
bool wavefront_reached(const Wavefront& wf, int32_t diag, int32_t anti_diag) {
    if (diag >= wf.diag_begin && diag < wf.diag_begin + wf.entries.size()) {
        return (wf.entries[diag - wf.diag_begin].M == anti_diag);
    }
    return false;
}

// always begins at the last wave front, creates the CIGAR in reverse order,
// appends new CIGAR operations to the CIGAR string that is passed in
template<typename WFVector>
void wavefront_traceback_internal(const std::string& seq1, const std::string& seq2,
                                  const WFScores& scores, const WFVector& wfs,
                                  int64_t& d, int64_t& lead_matches, WFMatrix_t& mat, int64_t& s,
                                  std::vector<CIGAROp>& cigar) {
    
    s = wfs.size() - 1;
    
    // TODO: traceback is O(N + M) even with suffix tree, but backpointers
    // could make it O(s). should i do that?
    
    uint32_t op_len = 0;
    while (true) {
        // we're in the match/mismatch matrix
        const auto& wf = wfs[s];
        if (mat == MAT_M) {
            int64_t a = wf.entries[d - wf.diag_begin].M;
            // TODO: i don't love this solution for the partial traceback problem...
            a -= 2 * lead_matches;
            lead_matches = 0;
            // extend through any matches
            int64_t i = (d + a) / 2;
            int64_t j = (a - d) / 2;
            while (i >= 0 && j >= 0 && seq1[i] == seq2[j]
                   && a != wf.entries[d - wf.diag_begin].I && a != wf.entries[d - wf.diag_begin].D) {
                op_len += 1;
                --i;
                --j;
                a -= 2;
                ++lead_matches;
            }
            //std::cerr << "matrix M, s " << s << ", match " << lead_matches << endl;
            if (s == 0) {
                // this condition handles the initial wf_extend, which was not preceded
                // by a wf_next, so we skip that part of the traceback
                break;
            }
            if (a == wf.entries[d - wf.diag_begin].I) {
                // this is where an insertion closed
                if (op_len != 0) {
                    cigar.emplace_back('M', op_len);
                }
                mat = MAT_I;
                op_len = 0;
                lead_matches = 0;
                continue;
            }
            else if (a == wf.entries[d - wf.diag_begin].D) {
                // this is where a deletion closed
                if (op_len != 0) {
                    cigar.emplace_back('M', op_len);
                }
                mat = MAT_D;
                op_len = 0;
                lead_matches = 0;
                continue;
            }
            else if (s >= scores.mismatch) {
                // this must a mismatch
                op_len += 1;
                s -= scores.mismatch;
                lead_matches = 0;
                continue;
            }
            // the traceback goes outside of the DP structure (should only happen
            // in the low-memory code path)
            break;
        }
        else if (mat == MAT_I) {
            // we're in the insertion matrix
            if (s >= scores.gap_extend) {
                const auto& wf_prev = wfs[s - scores.gap_extend];
                if (d - 1 >= wf_prev.diag_begin && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    if (wf.entries[d - wf.diag_begin].I == wf_prev.entries[d - 1 - wf_prev.diag_begin].I + 1) {
                        // an insert extended here
                        s -= scores.gap_extend;
                        op_len += 1;
                        d -= 1;
                        continue;
                    }
                }
            }
            if (s >= scores.gap_extend + scores.gap_open) {
                const auto& wf_prev = wfs[s - scores.gap_extend - scores.gap_open];
                if (d - 1 >= wf_prev.diag_begin  && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    if (wf.entries[d - wf.diag_begin].I == wf_prev.entries[d - 1 - wf_prev.diag_begin].M + 1) {
                        // an insert opened here
                        s -= scores.gap_extend + scores.gap_open;
                        op_len += 1;
                        cigar.emplace_back('I', op_len);
                        d -= 1;
                        mat = MAT_M;
                        op_len = 0;
                        continue;
                    }
                }
            }
            // the traceback goes outside of the DP structure (should only happen
            // in the low-memory code path)
            break;
        }
        else {
            // we're in the deletion matrix
            if (s >= scores.gap_extend) {
                const auto& wf_prev = wfs[s - scores.gap_extend];
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()) {
                    if (wf.entries[d - wf.diag_begin].D == wf_prev.entries[d + 1 - wf_prev.diag_begin].D + 1) {
                        // a deletion extended here
                        s -= scores.gap_extend;
                        op_len += 1;
                        d += 1;
                        continue;
                    }
                }
            }
            if (s >= scores.gap_extend + scores.gap_open) {
                const auto& wf_prev = wfs[s - scores.gap_extend - scores.gap_open];
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.entries.size()){
                    if (wf.entries[d - wf.diag_begin].D == wf_prev.entries[d + 1 - wf_prev.diag_begin].M + 1) {
                        // a deletion opened here
                        s -= scores.gap_extend + scores.gap_open;
                        op_len += 1;
                        d += 1;
                        cigar.emplace_back('D', op_len);
                        mat = MAT_M;
                        op_len = 0;
                        continue;
                    }
                }
            }
            // the traceback goes outside of the DP structure (should only happen
            // in the low-memory code path)
            break;
        }
    }
    
    // handle the final operation
    if (op_len != 0) {
        if (mat == MAT_M) {
            cigar.emplace_back('M', op_len);
        }
        else if (mat == MAT_D) {
            cigar.emplace_back('D', op_len);
        }
        else {
            cigar.emplace_back('I', op_len);
        }
    }
}

inline
std::vector<CIGAROp> wavefront_traceback(const std::string& seq1, const std::string& seq2,
                                         const WFScores& scores, const std::vector<Wavefront>& wfs) {
    
    // begin traceback at the end of the sequence
    int64_t d = seq1.size() - seq2.size();
    WFMatrix_t mat = MAT_M;
    int64_t dummy_s;
    int64_t lead_matches = 0;
    
    std::vector<CIGAROp> cigar;
    wavefront_traceback_internal(seq1, seq2, scores, wfs, d, lead_matches, mat, dummy_s, cigar);
    assert(dummy_s == 0 && d == 0 && mat == MAT_M);
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

namespace debug {

inline
void print_wf(const Wavefront& wf, int diag_begin, int diag_end) {
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.entries.size()) {
            int a = wf.entries[d - wf.diag_begin].M;
            if (a > -100) {
                std::cerr << a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.entries.size()) {
            int a = wf.entries[d - wf.diag_begin].I;
            if (a > -100) {
                std::cerr << a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.entries.size()) {
            int a = wf.entries[d - wf.diag_begin].D;
            if (a > -100) {
                std::cerr << a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
}

template<typename WFVector>
void print_wfs(const WFVector& wfs) {
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int s = 0; s < wfs.size(); ++s) {
        min_begin = std::min<int>(wfs[s].diag_begin, min_begin);
        max_end = std::max<int>(wfs[s].diag_begin + wfs[s].entries.size(), max_end);
    }
    std::cerr << "# print wfs #" << std::endl;
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int s = 0; s < wfs.size(); ++s) {
        std::cerr << "score " << s << std::endl;
        print_wf(wfs[s], min_begin, max_end);
    }
}

inline
void print_wfs_low_mem(const std::vector<std::pair<int32_t, Wavefront>>& wf_bank,
                       const std::deque<Wavefront>& wf_buffer, int s) {
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int i = 0; i < wf_bank.size(); ++i) {
        min_begin = std::min<int>(wf_bank[i].second.diag_begin, min_begin);
        max_end = std::max<int>(wf_bank[i].second.diag_begin + wf_bank[i].second.entries.size(), max_end);
    }
    for (int i = 0; i < wf_buffer.size(); ++i) {
        min_begin = std::min<int>(wf_buffer[i].diag_begin, min_begin);
        max_end = std::max<int>(wf_buffer[i].diag_begin + wf_buffer[i].entries.size(), max_end);
    }
    std::cerr << "### print wf bank ###" << std::endl;
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int t = 0, i = 0; !wf_bank.empty() && t <= wf_bank.back().first; ++t) {
        std::cerr << "score " << t << std::endl;
        while (t > wf_bank[i].first) {
            ++i;
        }
        if (t == wf_bank[i].first) {
            print_wf(wf_bank[i].second, min_begin, max_end);
        }
        else {
            for (int m = 0; m < 3; ++m) {
                for (int d = min_begin; d < max_end; ++d) {
                    std::cerr << "-\t";
                }
                std::cerr << std::endl;
            }
        }
    }
    std::cerr << "### print wf buffer ###" << std::endl;
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int i = 0; i < wf_buffer.size(); ++i) {
        std::cerr << "score " << s + i - wf_buffer.size() + 1 << std::endl;
        print_wf(wf_buffer[i], min_begin, max_end);
    }
}

inline
void print_wfs_tb(const std::deque<Wavefront>& traceback_block,
                  const std::vector<std::pair<int32_t, Wavefront>>& wf_bank,
                  int i) {
    
    std::cerr << "### print traceback block " << wf_bank[i].first << ":" << wf_bank[i].first + traceback_block.size() << " ###" << std::endl;
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int k = 0; k < traceback_block.size(); ++k) {
        min_begin = std::min<int>(min_begin, traceback_block[k].diag_begin);
        max_end = std::max<int>(max_end, traceback_block[k].diag_begin + traceback_block[k].entries.size());
    }
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int k = 0; k < traceback_block.size(); ++k) {
        std::cerr << "score " << wf_bank[i].first + k << std::endl;
        print_wf(traceback_block[k], min_begin, max_end);
    }
}

inline
void wavefront_viz(const std::string& seq1, const std::string& seq2,
                   const std::vector<Wavefront>& wfs) {
    
    std::vector<std::vector<int>> matrix(seq1.size() + 1, std::vector<int32_t>(seq2.size() + 1, -1));
    
    for (int s = 0; s < wfs.size(); ++s) {
        const auto& wf = wfs[s];
        for (int k = 0; k < wf.entries.size(); ++k) {
            int d = wf.diag_begin + k;
            int a = wf.entries[k].M;
            if (a < -2) {
                continue;
            }
            int i = (a + d) / 2;
            int j = (a - d) / 2;
            while (i >= 0 && j >= 0 && i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]
                   && a != wf.entries[k].I && a != wf.entries[k].D) {
                matrix[i + 1][j + 1] = s;
                --i;
                --j;
                a -= 2;
            }
            matrix[i + 1][j + 1] = s;
        }
    }
    
    std::cerr << "\t";
    for (auto c : seq2) {
        std::cerr << "\t" << c;
    }
    std::cerr << std::endl;
    for (int i = 0; i < seq1.size() + 1; ++i) {
        if (i != 0) {
            std::cerr << seq1[i - 1];
        }
        for (int j = 0; j < seq2.size() + 1; ++j) {
            std::cerr << '\t';
            if (matrix[i][j] >= 0) {
                std::cerr << matrix[i][j];
            }
            else {
                std::cerr << '.';
            }
        }
        std::cerr << std::endl;
    }
}

} // end namespace debug

inline
std::vector<CIGAROp> wavefront_traceback_low_mem(const std::string& seq1, const std::string& seq2,
                                                 const WFScores& scores, int32_t prune_diff,
                                                 const std::vector<std::pair<int32_t, Wavefront>>& wf_bank) {
    
    std::vector<CIGAROp> cigar;
    
    // a moving window over the DP that we will incrementally do traceback through
    std::deque<Wavefront> traceback_block;
    // a pointer to the first position in wf_bank that has been moved into
    // the traceback block
    int64_t i = wf_bank.size() - 1;
    traceback_block.emplace_front(std::move(wf_bank[i].second));
    
    while (i - 1 >= 0 && wf_bank[i - 1].first == wf_bank[i].first - 1) {
        // the scores are contiguous
        --i;
        traceback_block.emplace_front(std::move(wf_bank[i].second));
    }
    // begin traceback at the end of the sequence
    int64_t d = seq1.size() - seq2.size();
    WFMatrix_t mat = MAT_M;
    int64_t lead_matches = 0;
    // the index within the traceback block that traceback ends, initial value is arbitrary
    int64_t relative_s = 0;
//    std::cerr << "### begin traceback at diag " << d << " ###" << std::endl;
//    print_wfs_tb(traceback_block, wf_bank, i);
    
    // traceback as far as possible in this block
    wavefront_traceback_internal(seq1, seq2, scores, traceback_block,
                                 d, lead_matches, mat, relative_s, cigar);
    // remove anything we don't need anymore
    traceback_block.resize(relative_s + 1);
//    std::cerr << "conclude at score " << wf_bank[i].first + relative_s << ", diag " << d << " in matrix " << (mat == MAT_M ? "M" : (mat == MAT_I ? "I" : "D")) << std::endl;
//    std::cerr << "state of (reversed) cigar:" << std::endl;
//    for (auto cigar_op : cigar) {
//        std::cerr << cigar_op.len << cigar_op.op;
//    }
//    std::cerr << std::endl;
    
    while (i > 0) {
        // there's an earlier discontiguous stripe we need to connect to
        
        // figure out where it starts
        int64_t j = i - 1;
        while (j - 1 >= 0 && wf_bank[j - 1].first == wf_bank[j].first - 1) {
            --j;
        }
        
        // initialize the full block with stripe
        std::vector<Wavefront> fillin_wfs;
        fillin_wfs.reserve(wf_bank[i].first - wf_bank[j].first);
        for (int64_t k = j; k < i; ++k) {
            fillin_wfs.emplace_back(std::move(wf_bank[k].second));
        }
        
        // redo the DP between the stripes
        for (int32_t t = wf_bank[i - 1].first + 1, s = wf_bank[i].first; t < s; ++t) {
            fillin_wfs.emplace_back(wavefront_next(seq1, seq2, scores, fillin_wfs));
            wavefront_extend(seq1, seq2, fillin_wfs.back());
            if (prune_diff >= 0) {
                // prune lagging diagonals
                wavefront_prune(fillin_wfs.back(), prune_diff);
            }
        }
        
        // and move the redone wavefronts into the block buffer for traceback
        for (int64_t k = fillin_wfs.size() - 1; k >= 0; --k) {
            traceback_block.emplace_front(std::move(fillin_wfs[k]));
        }
        i = j;
        
//        print_wfs_tb(traceback_block, wf_bank, i);
        
        // traceback as far as possible in this block
        wavefront_traceback_internal(seq1, seq2, scores, traceback_block,
                                     d, lead_matches, mat, relative_s, cigar);
        // remove anything we don't need anymore
        traceback_block.resize(relative_s + 1);
        
//        std::cerr << "conclude at score " << wf_bank[j].first + relative_s << ", diag " << d << " in matrix " << (mat == MAT_M ? "M" : (mat == MAT_I ? "I" : "D")) << std::endl;
//        std::cerr << "state of (reversed) cigar:" << std::endl;
//        for (auto cigar_op : cigar) {
//            std::cerr << cigar_op.len << cigar_op.op;
//        }
//        std::cerr << std::endl;
    }
    
    // we sometimes need to coalesce CIGAR operations that occurred across
    // the boundaries of blocks
    size_t into = 0;
    for (size_t j = 1; j < cigar.size(); ++j) {
        if (cigar[j].op == cigar[into].op) {
            // merge
            cigar[into].len += cigar[j].len;
        }
        else {
            // don't merge and move the cursor
            ++into;
            if (j != into) {
                // we merged earlier, need to compact
                cigar[into] = cigar[j];
            }
        }
    }
    // remove the excess CIGAR operations
    cigar.resize(into + 1);
    // the CIGAR is built backwards
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

} // end namespace internal

inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align(const std::string& seq1, const std::string& seq2,
                const WFScores& scores, int32_t prune_diff) {
    
    if (seq1.size() + seq2.size() > std::numeric_limits<int32_t>::max()) {
        std::stringstream strm;
        strm << "error: WFA implementation can only align sequences of combined length < "  << std::numeric_limits<int32_t>::max() << '\n';
        throw std::runtime_error(strm.str());
    }
    
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size() - 2;
    
    // init wavefront to the upper left of both sequences' starts
    std::vector<internal::Wavefront> wfs(1, internal::Wavefront(0, 1));
    wfs.back().entries[0].M = -2;
    
    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wfs.front());
    while (!wavefront_reached(wfs.back(), final_diag, final_anti_diag)) {
        
        wfs.emplace_back(wavefront_next(seq1, seq2, scores, wfs));
        
        wavefront_extend(seq1, seq2, wfs.back());
        
        if (prune_diff >= 0) {
            // prune lagging diagonals
            wavefront_prune(wfs.back(), prune_diff);
        }
    }
    
#ifdef debug_viz
    wavefront_viz(seq1, seq2, wfs);
#endif
    
    return std::make_pair(wavefront_traceback(seq1, seq2, scores, wfs), int32_t(wfs.size() - 1));
}

inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem(const std::string& seq1, const std::string& seq2,
                        const WFScores& scores, int32_t prune_diff) {
    
    if (seq1.size() + seq2.size() > std::numeric_limits<int32_t>::max()) {
        std::stringstream strm;
        strm << "error: WFA implementation can only align sequences of combined length < "  << std::numeric_limits<int32_t>::max() << '\n';
        throw std::runtime_error(strm.str());
    }
    
    // epoch parameters used to sub-sample wavefronts to retain
    int32_t stripe_width = std::max(scores.mismatch,
                                    scores.gap_open + scores.gap_extend);
    int32_t epoch_len = 1;
    int32_t epoch_end = 1;
    int32_t sample_rate = 1;
    
    // use these to know when we've finished the alignment
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size() - 2;
    
    // init wavefront to the upper left of both sequences' starts
    std::deque<internal::Wavefront> wf_buffer(1, internal::Wavefront(0, 1));
    wf_buffer.back().entries[0].M = -2;
    int32_t s = 0;
    
    // the wavefronts we've decided to keep after they leave the buffer
    std::vector<std::pair<int32_t, internal::Wavefront>> wf_bank;
    
    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wf_buffer.front());
    while (!wavefront_reached(wf_buffer.back(), final_diag, final_anti_diag)) {
        
        
        // increase score and get the next wavefront
        ++s;
        //std::cerr << "iteration " << s << "\n";
        auto wf_next = wavefront_next(seq1, seq2, scores, wf_buffer);
        
        if (wf_buffer.size() >= stripe_width) {
            // the stripe of the wavefront that is falling out of the buffer
            int32_t stripe_num = (s - stripe_width) / stripe_width;
            
            // is this the first wavefront of a new epoch?
            if (stripe_num == epoch_end && stripe_num * stripe_width == s - stripe_width) {
                // we're entering a new epoch, adust the sub-sampling parameters
                // accordingly
                epoch_len *= 4;
                epoch_end += epoch_len;
                sample_rate *= 2;
            }
            //std::cerr << "ejecting wf " << s - stripe_width << ", stripe " << stripe_num << " (width " << stripe_width << "), epoch of size " << epoch_len << " ending at stripe " << epoch_end << " with sample rate " << sample_rate << ", stripe remainder = " << (stripe_num % sample_rate) << "\n";
            
            // bit-hacky sub for % that works because sample_rate is a power of 2
            if ((stripe_num & (sample_rate - 1)) == 0) {
                // we're ejecting a stripe that's being retained,
                wf_bank.emplace_back(s - stripe_width, std::move(wf_buffer.front()));
            }
            
            wf_buffer.pop_front();
        }
        wf_buffer.emplace_back(std::move(wf_next));
        
        // follow matches
        wavefront_extend(seq1, seq2, wf_buffer.back());
        
        if (prune_diff >= 0) {
            // prune lagging diagonals
            wavefront_prune(wf_buffer.back(), prune_diff);
        }
        
        //        print_wfs_low_mem(wf_bank, wf_buffer, s);
    }
    
    // move the final wavefronts from the buffer to the bank
    for (size_t i = 0; i < wf_buffer.size(); ++i) {
        wf_bank.emplace_back(s - wf_buffer.size() + i + 1, std::move(wf_buffer[i]));
    }
    
    return std::make_pair(wavefront_traceback_low_mem(seq1, seq2, scores, prune_diff, wf_bank), s);
}

}


#endif /* wfa_lm.hpp */