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
#include <sdsl/csa_bitcompressed.hpp>


namespace wfalm {

/*
 * Class that performs WFA with the standard, low-memory, or recusive variants.
 * This class differs from WFAligner in that it uses a longest-common-extension
 * algorithm based on a suffix tree to compute matches. In most cases, this is
 * slower than the default comparison-based match computation.
 */
template<int NumPW>
class WFAlignerST : public WFAligner<NumPW> {
    
public:
    
    /// Initialize with WFA-style score parameters. On opening an insertion or
    /// deletion, *both* the gap open and gap extend penalties are applied.
    /// Scores returned by alignment methods will also be WFA-style.
    WFAlignerST(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    
    /// Initialize with Smith-Waterman-Gotoh-style score parameters. On opening an
    /// insertion or deletion, *both* the gap open and gap extend penalties are applied.
    /// Scores returned by alignment methods will also be Smith-Waterman-Gotoh-style.
    WFAlignerST(uint32_t match, uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    
    /// Default constructor. Performs edit distance alignment.
    WFAlignerST();
    
    /*****************************
     *  Configurable parameters  *
     *****************************/
    
    /// If set to a number >= 0, prunes diagonals that fall this many anti-diagonals
    /// behind the furthest-reaching diagonal. This can lead to suboptimal alignments,
    /// but it also can increase speed and reduce memory use. If set to a number < 0,
    /// no pruning is performed.
    using WFAligner<NumPW>::lagging_diagonal_prune;
    
    /***********************
     *  Alignment methods  *
     ***********************/
    
    /// Globally align two sequences using the fastest algorithm that will remain
    /// constrained to a given maximum memory, using a suffix tree to compute matches
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///   target_mem   The target maximum memory use in bytes
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_adaptive(const char* seq1, size_t len1,
                             const char* seq2, size_t len2,
                             uint64_t max_mem) const;
    
    
    /// Globally align two sequences in O(s log s) memory using the recursive WFA algorithm,
    /// using a suffix tree to compute matches.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_recursive(const char* seq1, size_t len1,
                              const char* seq2, size_t len2) const;
    
    /// Globally align two sequences in O(s^3/2) memory using the low-memory WFA algorithm,
    /// using a suffix tree to compute matches.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_low_mem(const char* seq1, size_t len1,
                            const char* seq2, size_t len2) const;
    
    /// Globally align two sequences in O(s^2) memory using the standard WFA algorithm,
    /// using a suffix tree to compute matches.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align(const char* seq1, size_t len1,
                    const char* seq2, size_t len2) const;
    
    /// Locally align two sequences from a seed using the adaptive WFA as an
    /// alignment engine. *Can only called if WFAligner was initialized with
    /// Smith-Waterman-Gotoh scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   max_mem          The target maximum memory use in bytes
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_adaptive(const char* seq1, size_t len1,
                                   const char* seq2, size_t len2,
                                   uint64_t max_mem,
                                   size_t anchor_begin_1, size_t anchor_end_1,
                                   size_t anchor_begin_2, size_t anchor_end_2,
                                   bool anchor_is_match = true) const;

    /// Locally align two sequences from a seed using the recursive WFA as an
    /// alignment engine. *Can only called if WFAligner was initialized with
    /// Smith-Waterman-Gotoh scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_recursive(const char* seq1, size_t len1,
                                    const char* seq2, size_t len2,
                                    size_t anchor_begin_1, size_t anchor_end_1,
                                    size_t anchor_begin_2, size_t anchor_end_2,
                                    bool anchor_is_match = true) const;

    /// Locally align two sequences from a seed using the low memory WFA as an
    /// alignment engine. *Can only called if WFAligner was initialized with
    /// Smith-Waterman-Gotoh scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_low_mem(const char* seq1, size_t len1,
                                  const char* seq2, size_t len2,
                                  size_t anchor_begin_1, size_t anchor_end_1,
                                  size_t anchor_begin_2, size_t anchor_end_2,
                                  bool anchor_is_match = true) const;

    /// Locally align two sequences from a seed using the standard WFA as an
    /// alignment engine. *Can only called if WFAligner was initialized with
    /// Smith-Waterman-Gotoh scoring parameters in the constructor.*
    ///
    /// Args:
    ///   seq1             First sequence to be aligned (need not be null-terminated)
    ///   len1             Length of first sequence to be aligned
    ///   seq2             Second sequence to be aligned (need not be null-terminated)
    ///   len2             Length of second sequence to be aligned
    ///   anchor_begin_1   First index of the anchor on seq1
    ///   anchor_end_1     Past-the-last index of the anchor on seq1
    ///   anchor_begin_2   First index of the anchor on seq2
    ///   anchor_end_2     Past-the-last index of the anchor on seq2
    ///   anchor_is_match  If false, the anchor sequence will be aligned, otherwise
    ///                    it will be assumed to be a match
    ///
    /// Return value:
    ///   A tuple consisting of:
    ///     - CIGAR string for alignment
    ///     - the alignment score
    ///     - a pair of indexes indicating the interval of aligned sequence on seq1
    ///     - a pair of indexes indicating the interval of aligned sequence on seq2
    inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local(const char* seq1, size_t len1,
                          const char* seq2, size_t len2,
                          size_t anchor_begin_1, size_t anchor_end_1,
                          size_t anchor_begin_2, size_t anchor_end_2,
                          bool anchor_is_match = true) const;
    
};
    
    
    
    

    /*******************************************************************
     *******************************************************************
     ***             Only internal methods below here              ***
     *******************************************************************
     *******************************************************************/
    
    
    
    
    
    
    
    
    
    
    

typedef sdsl::cst_sada<sdsl::csa_bitcompressed<>,
                       sdsl::lcp_dac<1>,
                       sdsl::bp_support_gg<sdsl::nearest_neighbour_dictionary<1>,
                                           sdsl::rank_support_v<>,
                                           sdsl::select_support_mcl<>, 64>,
                       sdsl::rank_support_v<10, 2>,
                       sdsl::select_support_mcl<10, 2>> SuffixTree;

// a match computing function based on the suffix tree
struct STMatchFunc {
public:
    STMatchFunc() = default;
    ~STMatchFunc() = default;
    
protected:
    
    template<typename StringType>
    void init(const StringType& seq1, const StringType& seq2) {
        
        offset = seq1.size() + 1;
        
        // encoding:
        //  a/A  1
        //  c/C  2
        //  g/G  3
        //  t/T  4
        //  else 5
        //  $    6
        static const uint8_t separator = 6;
        static const uint8_t st_encode[256] {
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 0
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 32
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            
            uint8_t(5), uint8_t(1), uint8_t(5), uint8_t(2), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(3),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 64
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(4), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            uint8_t(5), uint8_t(1), uint8_t(5), uint8_t(2), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(3),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 96
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(4), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 128
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 160
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 192
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), //
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), // 224
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5),
            uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5), uint8_t(5)  //
        };
        
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
    
    inline size_t lce(size_t i, size_t j) const {
        // in sdsl, suffix array indexes are 0-based, but leaf ranks are 1-based
        return suffix_tree.depth(suffix_tree.lca(suffix_tree.select_leaf(suffix_tree.csa.isa[i] + 1),
                                                 suffix_tree.select_leaf(suffix_tree.csa.isa[offset + j] + 1)));
    }
    
    SuffixTree suffix_tree;
    size_t offset;
};

class FwdSTMatchFunc : STMatchFunc {
public:
    
    FwdSTMatchFunc(const StringView& seq1, const StringView& seq2) {
        init<StringView>(seq1, seq2);
    }
    
    inline size_t operator()(size_t i, size_t j) const {
        return lce(i, j);
    }
};

class RevSTMatchFunc : STMatchFunc {
public:
    
    RevSTMatchFunc(const RevStringView& seq1, const RevStringView& seq2) {
        init<RevStringView>(seq1, seq2);
    }
    
    inline size_t operator()(size_t i, size_t j) const {
        return lce(i, j);
    }
};

template<int NumPW>
WFAlignerST<NumPW>::WFAlignerST() : WFAligner<NumPW>() {}

template<int NumPW>
WFAlignerST<NumPW>::WFAlignerST(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend)
    : WFAligner<NumPW>(mismatch, gap_open, gap_extend)
{
    // only need to delegate to WFAligner
}

template<int NumPW>
WFAlignerST<NumPW>::WFAlignerST(uint32_t match, uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend)
    : WFAligner<NumPW>(match, mismatch, gap_open, gap_extend)
{
    // only need to delegate to WFAligner
}

template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAlignerST<NumPW>::wavefront_align(const char* seq1, size_t len1,
                             const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return WFAligner<NumPW>::template wavefront_dispatch<false, WFAligner<NumPW>::StdMem, FwdSTMatchFunc>(str1, str2,
                                                                      std::numeric_limits<uint64_t>::max());
}

template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAlignerST<NumPW>::wavefront_align_low_mem(const char* seq1, size_t len1,
                                   const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return WFAligner<NumPW>::template wavefront_dispatch<false, WFAligner<NumPW>::LowMem, FwdSTMatchFunc>(str1, str2,
                                                                      std::numeric_limits<uint64_t>::max());
}

template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAlignerST<NumPW>::wavefront_align_recursive(const char* seq1, size_t len1,
                                     const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return WFAligner<NumPW>::template wavefront_dispatch<false, WFAligner<NumPW>::Recursive, FwdSTMatchFunc>(str1, str2,
                                                                         std::numeric_limits<uint64_t>::max());
}

template<int NumPW>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAlignerST<NumPW>::wavefront_align_adaptive(const char* seq1, size_t len1,
                                      const char* seq2, size_t len2,
                                      uint64_t max_mem) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return WFAligner<NumPW>::template wavefront_dispatch<false, WFAligner<NumPW>::AdaptiveMem, FwdSTMatchFunc>(str1, str2, max_mem);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAlignerST<NumPW>::wavefront_align_local(const char* seq1, size_t len1,
                                 const char* seq2, size_t len2,
                                 size_t anchor_begin_1, size_t anchor_end_1,
                                 size_t anchor_begin_2, size_t anchor_end_2,
                                 bool anchor_is_match) const {

    // configure the alignment and match functions
    typedef typename WFAligner<NumPW>::template StandardSemilocalWFA<RevSTMatchFunc, RevStringView> PrefWFA;
    typedef typename WFAligner<NumPW>::template StandardGlobalWFA<FwdSTMatchFunc, StringView> AnchorWFA;
    typedef typename WFAligner<NumPW>::template StandardSemilocalWFA<FwdSTMatchFunc, StringView> SuffWFA;

    return WFAligner<NumPW>::template wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                                              std::numeric_limits<uint64_t>::max(),
                                                                                              anchor_begin_1, anchor_end_1,
                                                                                              anchor_begin_2, anchor_end_2,
                                                                                              anchor_is_match);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAlignerST<NumPW>::wavefront_align_local_low_mem(const char* seq1, size_t len1,
                                         const char* seq2, size_t len2,
                                         size_t anchor_begin_1, size_t anchor_end_1,
                                         size_t anchor_begin_2, size_t anchor_end_2,
                                         bool anchor_is_match) const {

    // configure the alignment and match functions
    typedef typename WFAligner<NumPW>::template LowMemSemilocalWFA<RevSTMatchFunc, RevStringView> PrefWFA;
    typedef typename WFAligner<NumPW>::template LowMemGlobalWFA<FwdSTMatchFunc, StringView> AnchorWFA;
    typedef typename WFAligner<NumPW>::template LowMemSemilocalWFA<FwdSTMatchFunc, StringView> SuffWFA;

    return WFAligner<NumPW>::template wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                                              std::numeric_limits<uint64_t>::max(),
                                                                                              anchor_begin_1, anchor_end_1,
                                                                                              anchor_begin_2, anchor_end_2,
                                                                                              anchor_is_match);
}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAlignerST<NumPW>::wavefront_align_local_recursive(const char* seq1, size_t len1,
                                             const char* seq2, size_t len2,
                                             size_t anchor_begin_1, size_t anchor_end_1,
                                             size_t anchor_begin_2, size_t anchor_end_2,
                                             bool anchor_is_match) const {

    // configure the alignment and match functions
    typedef typename WFAligner<NumPW>::template RecursiveSemilocalWFA<RevSTMatchFunc, RevStringView> PrefWFA;
    typedef typename WFAligner<NumPW>::template RecursiveGlobalWFA<FwdSTMatchFunc, StringView> AnchorWFA;
    typedef typename WFAligner<NumPW>::template RecursiveSemilocalWFA<FwdSTMatchFunc, StringView> SuffWFA;

    return WFAligner<NumPW>::template wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                                              std::numeric_limits<uint64_t>::max(),
                                                                                              anchor_begin_1, anchor_end_1,
                                                                                              anchor_begin_2, anchor_end_2,
                                                                                              anchor_is_match);

}

template<int NumPW>
inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAlignerST<NumPW>::wavefront_align_local_adaptive(const char* seq1, size_t len1,
                                             const char* seq2, size_t len2,
                                             uint64_t max_mem,
                                             size_t anchor_begin_1, size_t anchor_end_1,
                                             size_t anchor_begin_2, size_t anchor_end_2,
                                             bool anchor_is_match) const {

    // configure the alignment and match functions
    typedef typename WFAligner<NumPW>::template AdaptiveSemilocalWFA<RevSTMatchFunc, RevStringView> PrefWFA;
    typedef typename WFAligner<NumPW>::template AdaptiveGlobalWFA<FwdSTMatchFunc, StringView> AnchorWFA;
    typedef typename WFAligner<NumPW>::template AdaptiveSemilocalWFA<FwdSTMatchFunc, StringView> SuffWFA;
    
    return WFAligner<NumPW>::template wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2, max_mem,
                                                                                              anchor_begin_1, anchor_end_1,
                                                                                              anchor_begin_2, anchor_end_2,
                                                                                              anchor_is_match);

}

} // end namespace wfalm


#endif /* wfa_lm_st.hpp */
