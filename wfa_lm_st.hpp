//
//  wfa_lm_st.hpp
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

#ifndef wfa_lm_st_hpp_included
#define wfa_lm_st_hpp_included

#include "wfa_lm.hpp"

#include <sdsl/construct.hpp>
#include <sdsl/cst_sada.hpp>
#include <sdsl/cst_sct3.hpp>


namespace wfalm {

// TODO: choose better defaults for (inverse) suffix array sampling

// Configurable for a space/speed tradeoff
//typedef sdsl::cst_cst3<> SuffixTree;
typedef sdsl::cst_sada<> SuffixTree;

/// Align two sequences using the wavefront alignment algorithm, returns a CIGAR
/// string.
/// Optionally performs pruning (if prune_diff >= 0) of diagonals that are the
/// indicated difference behind the leading diagonal, measured in antidiagonals.
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem_st(const std::string& seq1, const std::string& seq2,
                           const WFScores& scores, int32_t prune_diff = -1);

/// Same as above except with Smith-Waterman-Gotoh style scoring parameter
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem_st(const std::string& seq1, const std::string& seq2,
                           const SWGScores& scores, int32_t prune_diff = -1);

/// Same semantics as above, but marginally faster and with a higher memory complexity
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_st(const std::string& seq1, const std::string& seq2,
                   const WFScores& scores, int32_t prune_diff = -1);

/// Same as above except with Smith-Waterman-Gotoh style scoring parameter
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_st(const std::string& seq1, const std::string& seq2,
                   const SWGScores& scores, int32_t prune_diff = -1);


/// Output includes two pairs of indexes, which correspond to the aligned intervals of
/// the first and second sequence respectively
inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_low_mem_st(const std::string& seq1, const std::string& seq2,
                                 size_t anchor_begin_1, size_t anchor_end_1,
                                 size_t anchor_begin_2, size_t anchor_end_2,
                                 const SWGScores& scores, bool anchor_is_match = true);

inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_st(const std::string& seq1, const std::string& seq2,
                         size_t anchor_begin_1, size_t anchor_end_1,
                         size_t anchor_begin_2, size_t anchor_end_2,
                         const SWGScores& scores, bool anchor_is_match = true);




/*******************************************************************
 *******************************************************************
 ***             Only internal functions below here              ***
 *******************************************************************
 *******************************************************************/


namespace internal {


// encoding:
//  a/A  1
//  c/C  2
//  g/G  3
//  t/T  4
//  else 5
//  sep  6
static constexpr uint8_t separator = 6;
static constexpr uint8_t st_encode[256] {
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 0
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 32
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
    
    uint8_t(5), uint8_t(1), uint8_t(5), uint8_t(2), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(3), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 64
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(4), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
    uint8_t(5), uint8_t(1), uint8_t(5), uint8_t(2), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(3), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 96
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(4), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
    
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 128
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 160
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
    
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),// 192
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),//
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),// 224
    uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5)//
};

template<typename StringType>
struct STMatchFunc {
    
    STMatchFunc(const StringType& seq1, const StringType& seq2) : offset(seq1.size() + 1) {
        
        // TODO: i wish this were possible with an implicit copy, but looks like probably not
        // since the sdsl::store_to_file that backends the in-memory construction is only
        // specialized for string and char*
        //  - but also this gets us around the StringView vs. string issues in local alignment...
        
        // encode the joined string and the sentinel in a 3-bit alphabet
        std::string combined_seq(seq1.size() + seq2.size() + 1, 0);
        size_t i;
        for (i = 0; i < seq1.size(); ++i) {
            combined_seq[i] = st_encode[seq1[i]];
        }
        combined_seq[i] = separator;
        ++i;
        for (size_t j = 0; j < seq2.size(); ++i, ++j) {
            combined_seq[i] = st_encode[seq2[j]];
        }
        
        // construct in memory with 1-byte words
        sdsl::construct_im(suffix_tree, combined_seq, 1);
    }
    
    inline size_t operator()(size_t i, size_t j) const {
        // in sdsl, suffix array indexes are 0-based, but leaf ranks are 1-based
        return suffix_tree.depth(suffix_tree.lca(suffix_tree.select_leaf(suffix_tree.csa.isa[i] + 1),
                                                 suffix_tree.select_leaf(suffix_tree.csa.isa[offset + j] + 1)));
    }
    
    SuffixTree suffix_tree;
    size_t offset;
};

} // end namespace internal

// TODO: would be nice to do this with less duplication from the other header...

inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_st(const std::string& seq1, const std::string& seq2,
                   const WFScores& scores, int32_t prune_diff) {
    internal::STMatchFunc<std::string> match_func(seq1, seq2);
    return internal::wavefront_align_core<false>(seq1, seq2, scores, prune_diff, 0, match_func);
}


inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem_st(const std::string& seq1, const std::string& seq2,
                           const WFScores& scores, int32_t prune_diff) {
    internal::STMatchFunc<std::string> match_func(seq1, seq2);
    return internal::wavefront_align_low_mem_core<false>(seq1, seq2, scores, prune_diff, 0, match_func);
}


inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_st(const std::string& seq1, const std::string& seq2,
                   const SWGScores& scores, int32_t prune_diff) {
    
    // make WFA params
    WFScores wf_scores;
    bool even;
    std::tie(wf_scores, even) = internal::convert_score_params(scores);
    
    // do the alignment
    internal::STMatchFunc<std::string> match_func(seq1, seq2);
    auto result = internal::wavefront_align_core<false>(seq1, seq2, wf_scores, prune_diff, 0, match_func);
    
    // convert the score back to SWG params
    result.second = internal::convert_score(scores, seq1.size(), seq2.size(), even, result.second);
    return result;
}

/// Same as above except with Smith-Waterman-Gotoh style scoring parameter
inline
std::pair<std::vector<CIGAROp>, int32_t>
wavefront_align_low_mem_st(const std::string& seq1, const std::string& seq2,
                           const SWGScores& scores, int32_t prune_diff) {
    
    // make WFA params
    WFScores wf_scores;
    bool even;
    std::tie(wf_scores, even) = internal::convert_score_params(scores);
    
    // do the alignment
    internal::STMatchFunc<std::string> match_func(seq1, seq2);
    auto result = internal::wavefront_align_low_mem_core<false>(seq1, seq2, wf_scores, prune_diff, 0, match_func);
    
    // convert the score back to SWG params
    result.second = internal::convert_score(scores, seq1.size(), seq2.size(), even, result.second);
    return result;
    
}

inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_st(const std::string& seq1, const std::string& seq2,
                         size_t anchor_begin_1, size_t anchor_end_1,
                         size_t anchor_begin_2, size_t anchor_end_2,
                         const SWGScores& scores, bool anchor_is_match) {
    
    // configure the alignment and match functions
    typedef internal::StandardSemilocalWFA<internal::STMatchFunc<internal::RevStringView>, internal::RevStringView> PrefWFA;
    typedef internal::StandardGlobalWFA<internal::STMatchFunc<internal::StringView>, internal::StringView> AnchorWFA;
    typedef internal::StandardSemilocalWFA<internal::STMatchFunc<internal::StringView>, internal::StringView> SuffWFA;
    
    return internal::wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, seq2,
                                                                             anchor_begin_1, anchor_end_1,
                                                                             anchor_begin_2, anchor_end_2,
                                                                             scores, anchor_is_match);
}

inline
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
wavefront_align_local_low_mem_st(const std::string& seq1, const std::string& seq2,
                                 size_t anchor_begin_1, size_t anchor_end_1,
                                 size_t anchor_begin_2, size_t anchor_end_2,
                                 const SWGScores& scores, bool anchor_is_match) {
    
    // configure the alignment and match functions
    typedef internal::LowMemSemilocalWFA<internal::STMatchFunc<internal::RevStringView>, internal::RevStringView> PrefWFA;
    typedef internal::LowMemGlobalWFA<internal::STMatchFunc<internal::StringView>, internal::StringView> AnchorWFA;
    typedef internal::LowMemSemilocalWFA<internal::STMatchFunc<internal::StringView>, internal::StringView> SuffWFA;
    
    return internal::wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, seq2,
                                                                             anchor_begin_1, anchor_end_1,
                                                                             anchor_begin_2, anchor_end_2,
                                                                             scores, anchor_is_match);
}

}


#endif /* wfa_lm_st.hpp */
