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
#include <tuple>

namespace wfalm {

/*
 * Score parameters for wavefront alignment. On opening a insertion or
 * deletion, *both* the gap open and gap extend penalties are applied.
 */
struct WFScores {
    WFScores(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    WFScores() = default;
    
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
    SWGScores() = default;
    
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

// TODO: revamp pruning, i don't like this interface

/// Align two sequences using the wavefront alignment algorithm, returns a CIGAR
/// string.
/// Optionally performs pruning (if prune_diff >= 0) of diagonals that are the
/// indicated difference behind the leading diagonal, measured in antidiagonals.
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem(const std::string& seq1, const std::string& seq2,
                        const WFScores& scores, int32_t prune_diff = -1);

/// Same as above except with Smith-Waterman-Gotoh style scoring parameter
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem(const std::string& seq1, const std::string& seq2,
                        const SWGScores& scores, int32_t prune_diff = -1);

/// Same semantics as above, but marginally faster and with a higher memory complexity
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align(const std::string& seq1, const std::string& seq2,
                const WFScores& scores, int32_t prune_diff = -1);

/// Same as above except with Smith-Waterman-Gotoh style scoring parameter
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align(const std::string& seq1, const std::string& seq2,
                const SWGScores& scores, int32_t prune_diff = -1);


/// Output includes two pairs of indexes, which correspond to the aligned intervals of
/// the first and second sequence respectively
inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_low_mem(const std::string& seq1, const std::string& seq2,
                              size_t anchor_begin_1, size_t anchor_end_1,
                              size_t anchor_begin_2, size_t anchor_end_2,
                              const SWGScores& scores, bool anchor_is_match = true);

inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local(const std::string& seq1, const std::string& seq2,
                      size_t anchor_begin_1, size_t anchor_end_1,
                      size_t anchor_begin_2, size_t anchor_end_2,
                      const SWGScores& scores, bool anchor_is_match = true);




/*******************************************************************
 *******************************************************************
 ***             Only internal functions below here              ***
 *******************************************************************
 *******************************************************************/


namespace internal {

/*
 * Adapter to avoid string copying for when aligning subintervals
 */
class StringView {
public:
    StringView(const std::string& seq, size_t i, size_t len)
    : seq(seq), i(i), len(len)
    {
        if (i > seq.size() || i + len > seq.size()) {
            throw std::runtime_error("StringView bounds do not fit in string\n");
        }
    }
    
    inline char operator[](size_t idx) const {
        return seq[i + idx];
    }
    
    inline size_t size() const {
        return len;
    }
    
private:
    
    const std::string& seq;
    size_t i;
    size_t len;
};

/*
 * Adapter to avoid string copying for when aligning subintervals
 * that also reverses the string
 */
class RevStringView {
public:
    // note: interval is provided in forward direction
    RevStringView(const std::string& seq, size_t i, size_t len)
    : seq(seq), last(i + len - 1), len(len)
    {
        if (i > seq.size() || last + 1 > seq.size()) {
            throw std::runtime_error("RevStringView bounds do not fit in string\n");
        }
    }
    
    inline char operator[](size_t idx) const {
        return seq[last - idx];
    }
    
    inline size_t size() const {
        return len;
    }
    
private:
    
    const std::string& seq;
    size_t last;
    size_t len;
};

/* For use in traceback */
enum WFMatrix_t {MAT_M = 0, MAT_I = 1, MAT_D = 2};

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

template<typename StringType, typename MatchFunc>
void wavefront_extend(const StringType& seq1, const StringType& seq2,
                      Wavefront& wf, const MatchFunc& match_func) {
    
    for (int64_t k = 0; k < wf.entries.size(); ++k) {
        int64_t diag = wf.diag_begin + k;
        int64_t anti_diag = wf.entries[k].M;
        int64_t i = (diag + anti_diag) / 2 + 1;
        int64_t j = (anti_diag - diag) / 2 + 1;
        if (i >= 0 && j >= 0) {
            wf.entries[k].M += 2 * match_func(i, j);
        }
    }
}

template<typename WFVector, typename StringType>
Wavefront wavefront_next(const StringType& seq1, const StringType& seq2,
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
    // but it should be only ~s * o extra work
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
    
    // some twiddly code so that we can call ctor only once and preserve RVO
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
    
    if (lo < hi) {
        
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

// creates the CIGAR in reverse order,
// appends new CIGAR operations to the CIGAR string that is passed in
template<typename WFVector, typename StringType>
void wavefront_traceback_internal(const StringType& seq1, const StringType& seq2,
                                  const WFScores& scores, const WFVector& wfs,
                                  int64_t& d, int64_t& lead_matches, WFMatrix_t& mat, int64_t& s,
                                  std::vector<CIGAROp>& cigar) {
        
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

template <typename StringType>
std::vector<CIGAROp> wavefront_traceback(const StringType& seq1, const StringType& seq2,
                                         const WFScores& scores, const std::vector<Wavefront>& wfs,
                                         int64_t s, int64_t d) {
    
    // always begin traceback in the match matrix
    WFMatrix_t mat = MAT_M;
    int64_t lead_matches = 0;
    
    std::vector<CIGAROp> cigar;
    wavefront_traceback_internal(seq1, seq2, scores, wfs, d, lead_matches, mat, s, cigar);
    assert(s == 0 && d == 0 && mat == MAT_M);
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

namespace debug {

template<typename StringType>
std::string str(const StringType& str) {
    std::stringstream strm;
    for (size_t i = 0; i < str.size(); ++i) {
        strm << str[i];
    }
    return strm.str();
}

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

template <typename StringType>
void wavefront_viz(const StringType& seq1, const StringType& seq2,
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

// merge all adjacent, equivalent operations
void coalesce_cigar(std::vector<CIGAROp>& cigar) {
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
}

/*
 * Allows the WF bank to be used in wf_next
 */
struct WFBankAdapter {
    WFBankAdapter(const std::vector<std::pair<int32_t, Wavefront>>& wf_bank) : wf_bank(wf_bank) {
        
    }
    
    inline const Wavefront& operator[](size_t i) const {
        return wf_bank[i].second;
    };
    
    inline size_t size() const {
        return wf_bank.size();
    }
    
    const std::vector<std::pair<int32_t, Wavefront>>& wf_bank;
};

template <bool Local, typename StringType, typename MatchFunc>
std::vector<CIGAROp> wavefront_traceback_low_mem(const StringType& seq1, const StringType& seq2,
                                                 const WFScores& scores, int32_t prune_diff,
                                                 std::vector<std::pair<int32_t, Wavefront>>& wf_bank,
                                                 int64_t s, int64_t d, const MatchFunc& match_func) {
    
    // the return value
    std::vector<CIGAROp> cigar;
    
    if (Local) {
        // if we overshot the opt traceback location, clear the way for it to be
        // in the last position of the bank
        while (wf_bank.back().first > s) {
            wf_bank.pop_back();
        }
        // and recompute out to the opt if necessary
        WFBankAdapter adapter(wf_bank);
        while (wf_bank.back().first < s) {
            wf_bank.emplace_back(wf_bank.back().first + 1,
                                 wavefront_next(seq1, seq2, scores, adapter));
            wavefront_extend(seq1, seq2, wf_bank.back().second, match_func);
            if (prune_diff >= 0) {
                // prune lagging diagonals
                wavefront_prune(wf_bank.back().second, prune_diff);
            }
        }
        // the traceback starting location should now be in the final wavefront
        // in the bank, regardless of global/local
    }
    
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
    // begin traceback in the match matrix
    WFMatrix_t mat = MAT_M;
    int64_t lead_matches = 0;
    // the index within the traceback block that traceback ends, initial value is arbitrary
    int64_t relative_s = traceback_block.size() - 1;;
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
        for (int32_t t = wf_bank[i - 1].first + 1, end = wf_bank[i].first; t < end; ++t) {
            fillin_wfs.emplace_back(wavefront_next(seq1, seq2, scores, fillin_wfs));
            wavefront_extend(seq1, seq2, fillin_wfs.back(), match_func);
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
        relative_s = traceback_block.size() - 1;
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
    coalesce_cigar(cigar);
    // the CIGAR is built backwards
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

template<typename StringType>
void find_local_opt(const StringType& seq1, const StringType& seq2,
                    int64_t s, const Wavefront& wf, int32_t match_score,
                    int64_t& opt, int64_t& opt_diag, int64_t& opt_s, size_t& max_s) {
    
    for (int64_t k = 0; k < wf.entries.size(); ++k) {
        // the score of the corresponding local alignment
        // note: have to add 2 because i decided to start the DP at -2, sigh...
        auto a = wf.entries[k].M + 2;
        int32_t local_s = match_score * a - s;
        if (local_s > opt && a >= std::abs(wf.diag_begin + k)) {
            // we found a new best that is inside the matrix
            opt = local_s;
            opt_diag = wf.diag_begin + k;
            opt_s = s;
            // update our bound on how far out we need to look
            max_s = s + match_score * (2 * std::min(seq1.size(), seq2.size()) - a);
        }
    }
}

template<bool Local, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_core(const StringType& seq1, const StringType& seq2,
                     const WFScores& scores, int32_t prune_diff, int32_t match_score,
                     const MatchFunc& match_func) {
    
    if (seq1.size() + seq2.size() > std::numeric_limits<int32_t>::max()) {
        std::stringstream strm;
        strm << "error: WFA implementation can only align sequences of combined length < "  << std::numeric_limits<int32_t>::max() << '\n';
        throw std::runtime_error(strm.str());
    }
    
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size() - 2;
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    size_t s = 0;
    
    // init wavefront to the upper left of both sequences' starts
    std::vector<Wavefront> wfs(1, Wavefront(0, 1));
    wfs.back().entries[0].M = -2;
    
    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wfs.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, s, wfs.front(), match_score,
                       opt, opt_diag, opt_s, max_s);
    }
    while (!wavefront_reached(wfs.back(), final_diag, final_anti_diag) &&
           s <= max_s) {
        
        ++s;
        wfs.emplace_back(wavefront_next(seq1, seq2, scores, wfs));
        
        wavefront_extend(seq1, seq2, wfs.back(), match_func);
        if (Local) {
            find_local_opt(seq1, seq2, s, wfs.back(), match_score,
                           opt, opt_diag, opt_s, max_s);
        }
        
        if (prune_diff >= 0) {
            // prune lagging diagonals
            wavefront_prune(wfs.back(), prune_diff);
        }
        
//        std::cerr << "WF state after extend at s " << s << std::endl;
//        debug::print_wfs(wfs);
    }
    
#ifdef debug_viz
    wavefront_viz(seq1, seq2, wfs);
#endif
    
    int64_t traceback_s, traceback_d;
    if (Local) {
        traceback_d = opt_diag;
        traceback_s = opt_s;
    }
    else {
        traceback_d = seq1.size() - seq2.size();
        traceback_s = wfs.size() - 1;
    }
    
//    debug::print_wfs(wfs);
    
//    std::cerr << "beginning traceback" << std::endl;
    
    return std::make_pair(wavefront_traceback(seq1, seq2, scores, wfs,
                                              traceback_s, traceback_d),
                          int32_t(traceback_s));
}


template<bool Local, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem_core(const StringType& seq1, const StringType& seq2,
                             const WFScores& scores, int32_t prune_diff, int32_t match_score,
                             const MatchFunc& match_func) {
    
    if (seq1.size() + seq2.size() > std::numeric_limits<int32_t>::max()) {
        std::stringstream strm;
        strm << "error: WFA implementation can only align sequences of combined length < "  << std::numeric_limits<int32_t>::max() << '\n';
        throw std::runtime_error(strm.str());
    }
    
    // epoch parameters used to sub-sample wavefronts to retain
    int64_t stripe_width = std::max(scores.mismatch,
                                    scores.gap_open + scores.gap_extend);
    int64_t epoch_len = 1;
    int64_t epoch_end = 1;
    int64_t sample_rate = 1;
    
    // use these to know when we've finished the alignment
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size() - 2;
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    
    // init wavefront to the upper left of both sequences' starts
    std::deque<Wavefront> wf_buffer(1, Wavefront(0, 1));
    wf_buffer.back().entries[0].M = -2;
    size_t s = 0;
    
    // the wavefronts we've decided to keep after they leave the buffer
    std::vector<std::pair<int32_t, Wavefront>> wf_bank;
    
    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wf_buffer.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, s, wf_buffer.front(), match_score,
                       opt, opt_diag, opt_s, max_s);
    }
    while (!wavefront_reached(wf_buffer.back(), final_diag, final_anti_diag) &&
           s <= max_s) {
        
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
        wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
        
        if (Local) {
            find_local_opt(seq1, seq2, s, wf_buffer.back(), match_score,
                           opt, opt_diag, opt_s, max_s);
        }
        
        if (prune_diff >= 0) {
            // prune lagging diagonals
            wavefront_prune(wf_buffer.back(), prune_diff);
        }
    }
    
    // move the final wavefronts from the buffer to the bank
    for (size_t i = 0; i < wf_buffer.size(); ++i) {
        wf_bank.emplace_back(s - wf_buffer.size() + i + 1, std::move(wf_buffer[i]));
    }
    
    int64_t traceback_s, traceback_d;
    if (Local) {
        traceback_d = opt_diag;
        traceback_s = opt_s;
        s = opt_s;
    }
    else {
        traceback_d = seq1.size() - seq2.size();
        traceback_s = s;
    }
    
    return std::make_pair(wavefront_traceback_low_mem<Local>(seq1, seq2, scores, prune_diff,
                                                             wf_bank, traceback_s, traceback_d,
                                                             match_func),
                          s);
}

inline std::pair<WFScores, bool> convert_score_params(const SWGScores& scores) {
    std::pair<WFScores, bool> result;
    // we can make scores a little smaller if the match is even
    if (scores.match % 2 == 0) {
        result.second = true;
        result.first.mismatch = scores.match + scores.mismatch;
        result.first.gap_open = scores.gap_open;
        result.first.gap_extend = scores.gap_extend + scores.match / 2;
    }
    else {
        result.second = false;
        result.first.mismatch = 2 * (scores.match + scores.mismatch);
        result.first.gap_open = 2 * scores.gap_open;
        result.first.gap_extend = 2 * scores.gap_extend + scores.match;
    }
    return result;
}

inline int32_t convert_score(const SWGScores& scores, size_t len1, size_t len2,
                             bool even, int32_t score) {
    if (even) {
        return scores.match * (len1 + len2) - score;
    }
    else {
        return (scores.match * (len1 + len2) - score) / 2;
    }
}

std::pair<size_t, size_t> cigar_base_length(const std::vector<CIGAROp>& cigar) {
    size_t len1 = 0, len2 = 0;
    for (const auto& c : cigar) {
        switch (c.op) {
            case 'M':
            case '=':
            case 'X':
                len1 += c.len;
                len2 += c.len;
                break;
            case 'I':
                len1 += c.len;
            case 'D':
                len2 += c.len;
                
            default:
                break;
        }
    }
    return std::make_pair(len1, len2);
}

template<typename PrefSemilocalWFA, typename AnchorGlobalWFA, typename SuffSemilocalWFA>
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_core(const std::string& seq1, const std::string& seq2,
                           size_t anchor_begin_1, size_t anchor_end_1,
                           size_t anchor_begin_2, size_t anchor_end_2,
                           const SWGScores& scores, bool anchor_is_match) {
    // make WFA params
    WFScores wf_scores;
    bool even;
    std::tie(wf_scores, even) = convert_score_params(scores);
    
    // tail align the prefix
    RevStringView pref1(seq1, 0, anchor_begin_1);
    RevStringView pref2(seq2, 0, anchor_begin_2);
    
    auto pref_result = PrefSemilocalWFA()(pref1, pref2, wf_scores, -1, scores.match);
    
    std::vector<CIGAROp> cigar = std::move(pref_result.first);
    std::reverse(cigar.begin(), cigar.end());
    auto aligned_pref_len = cigar_base_length(cigar);
    int32_t pref_score = convert_score(scores, aligned_pref_len.first, aligned_pref_len.second,
                                                 even, pref_result.second);
    
    int32_t anchor_score = 0;
    if (!anchor_is_match) {
        
        // align the anchor
        StringView anchor1(seq1, anchor_begin_1, anchor_end_1 - anchor_begin_1);
        StringView anchor2(seq2, anchor_begin_2, anchor_end_2 - anchor_begin_2);
        
        auto anchor_res = AnchorGlobalWFA()(anchor1, anchor2, wf_scores, -1, 0);
        
        // get the score and add to the CIGAR
        anchor_score = convert_score(scores, anchor_end_1 - anchor_begin_1,
                                               anchor_end_2 - anchor_begin_2, even, anchor_res.second);
        for (auto& op : anchor_res.first) {
            cigar.emplace_back(std::move(op));
        }
    }
    else {
        // infer the anchor match
        auto len = anchor_end_1 - anchor_begin_1;
        if (len != anchor_end_2 - anchor_begin_2) {
            throw std::runtime_error("error: match anchor for local alignment has mismatched lengths");
        }
        anchor_score = scores.match * len;
        cigar.emplace_back('M', len);
    }
    
    // tail align the suffix
    StringView suff1(seq1, anchor_end_1, seq1.size() - anchor_end_1);
    StringView suff2(seq2, anchor_end_2, seq2.size() - anchor_end_2);
    
    auto suff_result = SuffSemilocalWFA()(suff1, suff2, wf_scores, -1, scores.match);
    
    auto aligned_suff_len = cigar_base_length(suff_result.first);
    int32_t suff_score = convert_score(scores, aligned_suff_len.first, aligned_suff_len.second,
                                       even, suff_result.second);
    
    for (auto& op : suff_result.first) {
        cigar.emplace_back(std::move(op));
    }
    
    // CIGAR operations might not be merged over the boundaries
    coalesce_cigar(cigar);
    return std::make_tuple(std::move(cigar),
                           pref_score + anchor_score + suff_score,
                           std::make_pair(anchor_begin_1 - aligned_pref_len.first,
                                          anchor_end_1 + aligned_suff_len.first),
                           std::make_pair(anchor_begin_2 - aligned_pref_len.second,
                                          anchor_end_2 + aligned_suff_len.second));
}



template<typename StringType>
struct CompareMatchFunc {
    CompareMatchFunc(const StringType& seq1, const StringType& seq2)
    : seq1(seq1), seq2(seq2)
    {
        
    }
    
    inline size_t operator()(size_t i, size_t j) const {
        size_t init = i;
        while (i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]) {
            ++i;
            ++j;
        }
        return i - init;
    }
    
    const StringType& seq1;
    const StringType& seq2;
};


template<typename MatchFunc, typename StringType>
struct StandardGlobalWFA {
    StandardGlobalWFA() = default;
    
    inline std::pair<std::vector<CIGAROp>, int32_t>
    operator()(const StringType& seq1, const StringType& seq2,
               const WFScores& scores, int32_t prune_diff, int32_t match_score) const {
        MatchFunc match_func(seq1, seq2);
        return internal::wavefront_align_core<false>(seq1, seq2, scores, prune_diff,
                                                     match_score, match_func);
    }
};

template<typename MatchFunc, typename StringType>
struct StandardSemilocalWFA {
    StandardSemilocalWFA() = default;
    
    inline std::pair<std::vector<CIGAROp>, int32_t>
    operator()(const StringType& seq1, const StringType& seq2,
               const WFScores& scores, int32_t prune_diff, int32_t match_score) const {
        MatchFunc match_func(seq1, seq2);
        return internal::wavefront_align_core<true>(seq1, seq2, scores, prune_diff,
                                                    match_score, match_func);
    }
};

template<typename MatchFunc, typename StringType>
struct LowMemGlobalWFA {
    LowMemGlobalWFA() = default;
    
    inline std::pair<std::vector<CIGAROp>, int32_t>
    operator()(const StringType& seq1, const StringType& seq2,
               const WFScores& scores, int32_t prune_diff, int32_t match_score) const {
        MatchFunc match_func(seq1, seq2);
        return internal::wavefront_align_low_mem_core<false>(seq1, seq2, scores, prune_diff,
                                                             match_score, match_func);
    }
};

template<typename MatchFunc, typename StringType>
struct LowMemSemilocalWFA {
    LowMemSemilocalWFA() = default;
    
    inline std::pair<std::vector<CIGAROp>, int32_t>
    operator()(const StringType& seq1, const StringType& seq2,
               const WFScores& scores, int32_t prune_diff, int32_t match_score) const {
        MatchFunc match_func(seq1, seq2);
        return internal::wavefront_align_low_mem_core<true>(seq1, seq2, scores, prune_diff,
                                                            match_score, match_func);
    }
};

} // end namespace internal

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

inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align(const std::string& seq1, const std::string& seq2,
                const WFScores& scores, int32_t prune_diff) {
    internal::CompareMatchFunc<std::string> match_func(seq1, seq2);
    return internal::wavefront_align_core<false>(seq1, seq2, scores, prune_diff, 0, match_func);
}


inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem(const std::string& seq1, const std::string& seq2,
                        const WFScores& scores, int32_t prune_diff) {
    internal::CompareMatchFunc<std::string> match_func(seq1, seq2);
    return internal::wavefront_align_low_mem_core<false>(seq1, seq2, scores, prune_diff, 0, match_func);
}


inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align(const std::string& seq1, const std::string& seq2,
                const SWGScores& scores, int32_t prune_diff) {
    
    // make WFA params
    WFScores wf_scores;
    bool even;
    std::tie(wf_scores, even) = internal::convert_score_params(scores);
    
    // do the alignment
    internal::CompareMatchFunc<std::string> match_func(seq1, seq2);
    auto result = internal::wavefront_align_core<false>(seq1, seq2, wf_scores, prune_diff, 0, match_func);
    
    // convert the score back to SWG params
    result.second = internal::convert_score(scores, seq1.size(), seq2.size(), even, result.second);
    return result;
}

/// Same as above except with Smith-Waterman-Gotoh style scoring parameter
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem(const std::string& seq1, const std::string& seq2,
                        const SWGScores& scores, int32_t prune_diff) {
    
    // make WFA params
    WFScores wf_scores;
    bool even;
    std::tie(wf_scores, even) = internal::convert_score_params(scores);
    
    // do the alignment
    internal::CompareMatchFunc<std::string> match_func(seq1, seq2);
    auto result = internal::wavefront_align_low_mem_core<false>(seq1, seq2, wf_scores, prune_diff, 0, match_func);
    
    // convert the score back to SWG params
    result.second = internal::convert_score(scores, seq1.size(), seq2.size(), even, result.second);
    return result;
    
}

inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local(const std::string& seq1, const std::string& seq2,
                      size_t anchor_begin_1, size_t anchor_end_1,
                      size_t anchor_begin_2, size_t anchor_end_2,
                      const SWGScores& scores, bool anchor_is_match) {
    
    // configure the alignment and match functions
    typedef internal::StandardSemilocalWFA<internal::CompareMatchFunc<internal::RevStringView>, internal::RevStringView> PrefWFA;
    typedef internal::StandardGlobalWFA<internal::CompareMatchFunc<internal::StringView>, internal::StringView> AnchorWFA;
    typedef internal::StandardSemilocalWFA<internal::CompareMatchFunc<internal::StringView>, internal::StringView> SuffWFA;
    
    return internal::wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, seq2,
                                                                             anchor_begin_1, anchor_end_1,
                                                                             anchor_begin_2, anchor_end_2,
                                                                             scores, anchor_is_match);
}

inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_low_mem(const std::string& seq1, const std::string& seq2,
                              size_t anchor_begin_1, size_t anchor_end_1,
                              size_t anchor_begin_2, size_t anchor_end_2,
                              const SWGScores& scores, bool anchor_is_match) {
    
    // configure the alignment and match functions
    typedef internal::LowMemSemilocalWFA<internal::CompareMatchFunc<internal::RevStringView>, internal::RevStringView> PrefWFA;
    typedef internal::LowMemGlobalWFA<internal::CompareMatchFunc<internal::StringView>, internal::StringView> AnchorWFA;
    typedef internal::LowMemSemilocalWFA<internal::CompareMatchFunc<internal::StringView>, internal::StringView> SuffWFA;
    
    return internal::wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, seq2,
                                                                             anchor_begin_1, anchor_end_1,
                                                                             anchor_begin_2, anchor_end_2,
                                                                             scores, anchor_is_match);
}

}


#endif /* wfa_lm.hpp */
