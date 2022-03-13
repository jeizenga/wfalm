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
 * Operation in a CIGAR string (as defined in SAM format specification)
 */
struct CIGAROp {
    CIGAROp(char op, uint32_t len) : op(op), len(len) {}
    CIGAROp() = default;
    
    char op;
    uint32_t len;
};

/*
 * Class that performs WFA with the standard, low-memory, or recusive variants.
 * The aligner objects are relatively lightweight and can be constructed for
 * single-use with limited overhead.
 */
class WFAligner {
    
public:
    
    /// Initialize with WFA-style score parameters. On opening an insertion or
    /// deletion, *both* the gap open and gap extend penalties are applied.
    /// Scores returned by alignment methods will also be WFA-style.
    WFAligner(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    
    /// Initialize with Smith-Waterman-Gotoh-style score parameters. On opening an
    /// insertion or deletion, *both* the gap open and gap extend penalties are applied.
    /// Scores returned by alignment methods will also be Smith-Waterman-Gotoh-style.
    WFAligner(uint32_t match, uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend);
    
    /// Default constructor. Performs edit distance alignment.
    WFAligner();

    /*****************************
     *  Configurable parameters  *
     *****************************/
    
    /// If set to a number >= 0, prunes diagonals that fall this many anti-diagonals
    /// behind the furthest-reaching diagonal. This can lead to suboptimal alignments,
    /// but it also can increase speed and reduce memory use. If set to a number < 0,
    /// no pruning is performed.
    int32_t lagging_diagonal_prune = -1;
    
    /***********************
     *  Alignment methods  *
     ***********************/
    
    /// Globally align two sequences using the fastest algorithm that will remain
    /// constrained to a given maximum memory.
    ///
    /// Args:
    ///   seq1         First sequence to be aligned (need not be null-terminated)
    ///   len1         Length of first sequence to be aligned
    ///   seq2         Second sequence to be aligned (need not be null-terminated)
    ///   len2         Length of second sequence to be aligned
    ///   max_mem      The target maximum memory use in bytes
    ///
    /// Return value:
    ///   Pair consisting of CIGAR string for alignment and the alignment score.
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_adaptive(const char* seq1, size_t len1,
                             const char* seq2, size_t len2,
                             uint64_t max_mem) const;
    
    /// Globally align two sequences in O(s log s) memory and O(sN log s) time using the recursive WFA algorithm.
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
    
    /// Globally align two sequences in O(s^3/2) memory and O(sN) time using the low-memory WFA algorithm.
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
    
    /// Globally align two sequences in O(s^2) memory and O(sN) time using the standard WFA algorithm.
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
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
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
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
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
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
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
    /// alignment engine.
    /// *Can only called if WFAligner was initialized with Smith-Waterman-Gotoh
    /// scoring parameters in the constructor.*
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
    
    
    
    
    
    
    
    
    /*******************************************************************
     *******************************************************************
     ***             Only internal functions below here              ***
     *******************************************************************
     *******************************************************************/
    
    
    
    
    
    
    
    
    
protected:
    
    enum WFMatrix_t {MAT_M = 0, MAT_I = 1, MAT_D = 2};

    static const int StdMem = 0;
    static const int LowMem = 1;
    static const int Recursive = 2;
    static const int AdaptiveMem = 3;
    
    /*
     * Adapter to avoid string copying for when aligning subintervals
     */
    class StringView {
    public:
        StringView(const char* seq, size_t i, size_t len) : seq(seq + i), len(len)
        {
            
        }
        
        inline char operator[](size_t idx) const {
            return seq[idx];
        }
        
        inline size_t size() const {
            return len;
        }
        
    private:
        
        const char* seq;
        size_t len;
    };
    
    /*
     * Adapter to avoid string copying for when aligning subintervals
     * that also reverses the string
     */
    class RevStringView {
    public:
        // note: interval is provided in forward direction
        RevStringView(const char* seq, size_t i, size_t len) : seql(seq + i + len - 1), len(len)
        {
            
        }
        
        inline char operator[](size_t idx) const {
            return *(seql - idx);
        }
        
        inline size_t size() const {
            return len;
        }
        
    private:
        // last char in the seq
        const char* seql;
        size_t len;
    };
    
    // TODO: use posix_memalign?
    template<typename IntType>
    struct Wavefront {
        
    private:
        
        inline void init_arrays() {
            alloced = (IntType*) malloc(3 * len * sizeof(IntType));
            M = alloced;
            I = M + len;
            D = I + len;
        }
    public:
        
        Wavefront() {}
        Wavefront(int32_t diag_begin, int32_t diag_end)
            : diag_begin(diag_begin), len(diag_end - diag_begin)
        {
            init_arrays();
            for (int32_t i = 0, n = 3 * len; i < n; ++i) {
                alloced[i] = std::numeric_limits<IntType>::min();
            }
        }
        
        inline Wavefront& operator=(const Wavefront& other) noexcept {
            if (this != &other) {
                diag_begin = other.diag_begin;
                len = other.len;
                init_arrays();
                for (int32_t i = 0; i < len; ++i) {
                    M[i] = other.M[i];
                    I[i] = other.I[i];
                    D[i] = other.D[i];
                }
            }
            return *this;
        }
        
        inline Wavefront& operator=(Wavefront&& other) noexcept {
            if (this != &other) {
                diag_begin = other.diag_begin;
                len = other.len;
                if (alloced) {
                    free(alloced);
                }
                alloced = other.alloced;
                M = other.M;
                I = other.I;
                D = other.D;
                
                // TODO: i could probably get away with only resetting the
                // alloc'ed pointer...
                other.alloced = other.M = other.I = other.D = nullptr;
                other.len = 0;
            }
            return *this;
        }
        
        Wavefront(const Wavefront& other) noexcept {
            *this = other;
        }
        
        Wavefront(Wavefront&& other) noexcept {
            *this = std::move(other);
        }
        
        ~Wavefront() {
            free(alloced);
        }
        
        inline void pop_front() {
            ++M;
            ++I;
            ++D;
            ++diag_begin;
            --len;
        }
        
        inline void pop_back() {
            --len;
        }
        
        inline int64_t size() const {
            return len;
        }
        
        inline bool empty() const {
            return len == 0;
        }
        
        inline uint64_t bytes() {
            return sizeof(IntType) * 3 * len;
        }
        
        IntType* alloced = nullptr;
        IntType* M = nullptr;
        IntType* I = nullptr;
        IntType* D = nullptr;
        
        int32_t diag_begin = 0;
        int32_t len = 0;
    };
    
    
    /*
     * Allows the WF bank to be used in wf_next
     */
    template<typename IntType>
    struct WFBankAdapter {
        WFBankAdapter(const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank) : wf_bank(wf_bank) {
            
        }
        
        inline const Wavefront<IntType>& operator[](size_t i) const {
            return wf_bank[i].second;
        };
        
        inline size_t size() const {
            return wf_bank.size();
        }
        
        const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank;
    };
    
    // a match-computing function based on direct character comparison
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
    
    /*
     * wrappers to perform WFA alignment algorithms for dependency injection
     */
    template<typename MatchFunc, typename StringType>
    struct StandardGlobalWFA {
        StandardGlobalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner::StdMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                std::numeric_limits<uint64_t>::max());
        }
        const WFAligner* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct StandardSemilocalWFA {
        StandardSemilocalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner::StdMem, MatchFunc, StringType>(seq1, seq2,
                                                                                               std::numeric_limits<uint64_t>::max());
        }
        const WFAligner* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct LowMemGlobalWFA {
        LowMemGlobalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner::LowMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                std::numeric_limits<uint64_t>::max());
        }
        const WFAligner* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct LowMemSemilocalWFA {
        LowMemSemilocalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner::LowMem, MatchFunc, StringType>(seq1, seq2,
                                                                                               std::numeric_limits<uint64_t>::max());
        }
        const WFAligner* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct RecursiveGlobalWFA {
        RecursiveGlobalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner::Recursive, MatchFunc, StringType>(seq1, seq2,
                                                                                                   std::numeric_limits<uint64_t>::max());
        }
        const WFAligner* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct RecursiveSemilocalWFA {
        RecursiveSemilocalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner::Recursive, MatchFunc, StringType>(seq1, seq2,
                                                                                                  std::numeric_limits<uint64_t>::max());
        }
        const WFAligner* aligner = nullptr;
    };
    template<typename MatchFunc, typename StringType>
    struct AdaptiveGlobalWFA {
        AdaptiveGlobalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner), max_mem(max_mem) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<false, WFAligner::AdaptiveMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                     max_mem);
        }
        const WFAligner* aligner = nullptr;
        uint64_t max_mem = 0;
    };
    template<typename MatchFunc, typename StringType>
    struct AdaptiveSemilocalWFA {
        AdaptiveSemilocalWFA(const WFAligner* aligner, uint64_t max_mem) : aligner(aligner), max_mem(max_mem) {}
        
        inline std::pair<std::vector<CIGAROp>, int32_t>
        operator()(const StringType& seq1, const StringType& seq2) const {
            return aligner->wavefront_dispatch<true, WFAligner::AdaptiveMem, MatchFunc, StringType>(seq1, seq2,
                                                                                                    max_mem);
        }
        const WFAligner* aligner = nullptr;
        uint64_t max_mem = 0;
    };
    
    // give the wrapper access to the dispatch function
    friend class StandardGlobalWFA<CompareMatchFunc<StringView>, StringView>;
    friend class StandardSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView>;
    friend class LowMemGlobalWFA<CompareMatchFunc<StringView>, StringView>;
    friend class LowMemSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView>;
    friend class RecursiveGlobalWFA<CompareMatchFunc<StringView>, StringView>;
    friend class RecursiveSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView>;
    friend class AdaptiveGlobalWFA<CompareMatchFunc<StringView>, StringView>;
    friend class AdaptiveSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView>;
    
    // central routine that back-ends the user-facing interface
    template<bool Local, int Mem, typename MatchFunc, typename StringType>
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_dispatch(const StringType& seq1, const StringType& seq2, uint64_t max_mem) const;
    
    // greedy take matches in WFA
    template<typename StringType, typename MatchFunc, typename IntType>
    inline void wavefront_extend(const StringType& seq1, const StringType& seq2,
                                 Wavefront<IntType>& wf, const MatchFunc& match_func) const;
    
    // increase score and initialize row in WFA
    template<typename IntType, typename WFVector, typename StringType>
    inline Wavefront<IntType> wavefront_next(const StringType& seq1, const StringType& seq2,
                                             const WFVector& wfs) const;
    
    // prune lagging diagonals in WFA
    template<bool Local, typename IntType>
    inline void wavefront_prune(Wavefront<IntType>& wf) const;
    
    // check if reached complete alignment
    template<typename IntType>
    inline bool wavefront_reached(const Wavefront<IntType>& wf, int32_t diag, int32_t anti_diag) const;
    
    // traceback through a complete DP structure
    template <typename StringType, typename IntType>
    inline std::vector<CIGAROp> wavefront_traceback(const StringType& seq1, const StringType& seq2,
                                                    const std::vector<Wavefront<IntType>>& wfs,
                                                    int64_t s, int64_t d) const;
    
    // traceback as far as possible through a possibly incomplete DP structure
    // creates the CIGAR in reverse order,
    // appends new CIGAR operations to the CIGAR string that is passed in
    template<typename WFVector, typename StringType>
    void wavefront_traceback_internal(const StringType& seq1, const StringType& seq2,
                                      const WFVector& wfs, int64_t& d, int64_t& lead_matches,
                                      WFMatrix_t& mat, int64_t& s, std::vector<CIGAROp>& cigar) const;
    
    // traceback through the low memory algorithms DP structure, with recomputation
    template <bool Local, typename StringType, typename MatchFunc, typename IntType>
    std::vector<CIGAROp> wavefront_traceback_low_mem(const StringType& seq1, const StringType& seq2,
                                                     std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                                     int64_t s, int64_t d, const MatchFunc& match_func) const;
    
    // do classic WFA
    template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_core(const StringType& seq1, const StringType& seq2,
                         const MatchFunc& match_func, uint64_t max_mem) const;
    
    // do low memory WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_low_mem_core(const StringType& seq1, const StringType& seq2,
                                 const MatchFunc& match_func, uint64_t max_mim) const;
    
    // internal routine for low memory WFA after initialization or fall-back
    template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_low_mem_core_internal(const StringType& seq1, const StringType& seq2,
                                          const MatchFunc& match_func, uint64_t max_mem,
                                          int64_t epoch_len, int64_t epoch_end, int64_t sample_rate,
                                          int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                          std::deque<Wavefront<IntType>>& wf_buffer, size_t s,
                                          std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                          uint64_t curr_mem_buffer, uint64_t curr_mem_bank,
                                          uint64_t curr_mem_block, uint64_t max_mem_block) const;
    
    // do recursive WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_recursive_core(const StringType& seq1, const StringType& seq2,
                                   const MatchFunc& match_func) const;
    
    // internal routine for recursive WFA after initialization or fall-back
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    inline std::pair<std::vector<CIGAROp>, int32_t>
    wavefront_align_recursive_core_internal(const StringType& seq1, const StringType& seq2,
                                            const MatchFunc& match_func,
                                            int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                            std::vector<Wavefront<IntType>>& opt_stripe,
                                            std::vector<Wavefront<IntType>>& first_stripe,
                                            std::deque<Wavefront<IntType>>& wf_buffer, size_t s) const;
    
    // recursive calls of recursive WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc, typename WFVector1, typename WFVector2>
    void wavefront_align_recursive_internal(const StringType& seq1, const StringType& seq2,
                                            int64_t lower_stripe_num, WFVector1& lower_stripe,
                                            int64_t upper_stripe_num, WFVector2& upper_stripe,
                                            int64_t& traceback_s, int64_t& traceback_d,
                                            WFMatrix_t& traceback_mat, int64_t& lead_matches,
                                            std::vector<CIGAROp>& cigar,
                                            const MatchFunc& match_func) const;
    
    // give up on a standard WFA and switch to doing low memory WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    fall_back_to_low_mem(const StringType& seq1, const StringType& seq2,
                         const MatchFunc& match_func, uint64_t max_mem,
                         int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                         std::vector<Wavefront<IntType>>& wfs) const;
    
    // give up on low memory WFA and switch to doing recursive WFA
    template<bool Local, typename IntType, typename StringType, typename MatchFunc>
    std::pair<std::vector<CIGAROp>, int32_t>
    fall_back_to_recursive(const StringType& seq1, const StringType& seq2,
                           const MatchFunc& match_func,
                           int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                           std::deque<Wavefront<IntType>>& wf_buffer, size_t s,
                           std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank) const;
    
    // check score that are SWG-locally optimal
    template<typename StringType, typename IntType>
    inline void find_local_opt(const StringType& seq1, const StringType& seq2,
                               int64_t s, const Wavefront<IntType>& wf,
                               int64_t& opt, int64_t& opt_diag, int64_t& opt_s, size_t& max_s) const;
    
    // do anchored local alignment
    template<typename PrefSemilocalWFA, typename AnchorGlobalWFA, typename SuffSemilocalWFA>
    std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
    wavefront_align_local_core(const char* seq1, size_t len1,
                               const char* seq2, size_t len2,
                               uint64_t max_mem,
                               size_t anchor_begin_1, size_t anchor_end_1,
                               size_t anchor_begin_2, size_t anchor_end_2,
                               bool anchor_is_match) const;
    
    // greatest commmon denominator
    uint32_t gcd(uint32_t a, uint32_t b) const;
    
    // convert a WFA score into a SWG score
    inline int32_t convert_score(size_t len1, size_t len2, int32_t score) const;
    
    // merge adjacent CIGAR ops that are the same type
    inline void coalesce_cigar(std::vector<CIGAROp>& cigar) const;
    
    // get the length of the aligned sequence on both sequences
    inline std::pair<size_t, size_t> cigar_base_length(const std::vector<CIGAROp>& cigar) const;
    
    // WFA-style parameters
    uint32_t mismatch = 1;
    uint32_t gap_open = 0;
    uint32_t gap_extend = 1;
    
    // scores are reduced by a constant factor if they have a non-trivial common
    // denominator to reduce run time
    uint32_t factor = 1;
    
    // SWG match score, or 0 if not using SWG scoring
    uint32_t match = 0;
    
    // how many subsequent wavefronts are required to do DP with these scoring parameters
    int64_t stripe_width = 1;
    
private:
    
    /* debug functions*/
    template<typename IntType>
    inline void print_wf(const Wavefront<IntType>& wf, int diag_begin, int diag_end) const;
    template<typename WFVector>
    void print_wfs(const WFVector& wfs) const;
    template<typename IntType>
    inline void print_wfs_low_mem(const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                  const std::deque<Wavefront<IntType>>& wf_buffer, int s) const;
    template<typename IntType>
    inline void print_wfs_tb(const std::deque<Wavefront<IntType>>& traceback_block,
                             const std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                             int i) const;
    template <typename StringType, typename WFVector>
    void wavefront_viz(const StringType& seq1, const StringType& seq2,
                       const WFVector& wfs) const;
};

// default
inline WFAligner::WFAligner() {}

inline WFAligner::WFAligner(uint32_t _mismatch, uint32_t _gap_open, uint32_t _gap_extend) {
    // reduce by common factor to accelerate the alignment
    if (_mismatch == 0 || _gap_extend == 0) {
        std::stringstream strm;
        strm << "error:[WFAligner] mismatch and extend penalties must both be non-zero" << std::endl;
        throw std::runtime_error(strm.str());
    }
    factor = gcd(_mismatch, _gap_extend);
    if (_gap_open != 0) {
        factor = gcd(factor, _gap_open);
    }
    mismatch = _mismatch / factor;
    gap_open = _gap_open / factor;
    gap_extend = _gap_extend / factor;
    
    stripe_width = std::max(gap_open + gap_extend, mismatch);
};

inline WFAligner::WFAligner(uint32_t _match, uint32_t _mismatch, uint32_t _gap_open, uint32_t _gap_extend)
    : WFAligner(2 * (_match + _mismatch), 2 * _gap_open, 2 * _gap_extend + _match)
{
    if (_match == 0) {
        std::stringstream strm;
        strm << "error:[WFAligner] match score must be an integer > 0 for Smith-Waterman-Gotoh-style score parameters: " << match << "\n";
        throw std::runtime_error(strm.str());
    }
    match = _match;
}


inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align(const char* seq1, size_t len1,
                           const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, StdMem, CompareMatchFunc<StringView>>(str1, str2,
                                                                           std::numeric_limits<uint64_t>::max());
}


inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_low_mem(const char* seq1, size_t len1,
                                   const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, LowMem, CompareMatchFunc<StringView>>(str1, str2,
                                                                           std::numeric_limits<uint64_t>::max());
}


inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_recursive(const char* seq1, size_t len1,
                                     const char* seq2, size_t len2) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, Recursive, CompareMatchFunc<StringView>>(str1, str2,
                                                                              std::numeric_limits<uint64_t>::max());
}

inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_adaptive(const char* seq1, size_t len1,
                                    const char* seq2, size_t len2,
                                    uint64_t max_mem) const {
    StringView str1(seq1, 0, len1);
    StringView str2(seq2, 0, len2);
    return wavefront_dispatch<false, AdaptiveMem, CompareMatchFunc<StringView>>(str1, str2, max_mem);
}

inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner::wavefront_align_local(const char* seq1, size_t len1,
                                 const char* seq2, size_t len2,
                                 size_t anchor_begin_1, size_t anchor_end_1,
                                 size_t anchor_begin_2, size_t anchor_end_2,
                                 bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef StandardSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView> PrefWFA;
    typedef StandardGlobalWFA<CompareMatchFunc<StringView>, StringView> AnchorWFA;
    typedef StandardSemilocalWFA<CompareMatchFunc<StringView>, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                   std::numeric_limits<uint64_t>::max(),
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner::wavefront_align_local_low_mem(const char* seq1, size_t len1,
                                         const char* seq2, size_t len2,
                                         size_t anchor_begin_1, size_t anchor_end_1,
                                         size_t anchor_begin_2, size_t anchor_end_2,
                                         bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef LowMemSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView> PrefWFA;
    typedef LowMemGlobalWFA<CompareMatchFunc<StringView>, StringView> AnchorWFA;
    typedef LowMemSemilocalWFA<CompareMatchFunc<StringView>, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                   std::numeric_limits<uint64_t>::max(),
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner::wavefront_align_local_recursive(const char* seq1, size_t len1,
                                           const char* seq2, size_t len2,
                                           size_t anchor_begin_1, size_t anchor_end_1,
                                           size_t anchor_begin_2, size_t anchor_end_2,
                                           bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef RecursiveSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView> PrefWFA;
    typedef RecursiveGlobalWFA<CompareMatchFunc<StringView>, StringView> AnchorWFA;
    typedef RecursiveSemilocalWFA<CompareMatchFunc<StringView>, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2,
                                                                   std::numeric_limits<uint64_t>::max(),
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

inline std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner::wavefront_align_local_adaptive(const char* seq1, size_t len1,
                                          const char* seq2, size_t len2,
                                          uint64_t max_mem,
                                          size_t anchor_begin_1, size_t anchor_end_1,
                                          size_t anchor_begin_2, size_t anchor_end_2,
                                          bool anchor_is_match) const {
    
    // configure the alignment and match functions
    typedef AdaptiveSemilocalWFA<CompareMatchFunc<RevStringView>, RevStringView> PrefWFA;
    typedef AdaptiveGlobalWFA<CompareMatchFunc<StringView>, StringView> AnchorWFA;
    typedef AdaptiveSemilocalWFA<CompareMatchFunc<StringView>, StringView> SuffWFA;
    
    return wavefront_align_local_core<PrefWFA, AnchorWFA, SuffWFA>(seq1, len1, seq2, len2, max_mem,
                                                                   anchor_begin_1, anchor_end_1,
                                                                   anchor_begin_2, anchor_end_2,
                                                                   anchor_is_match);
}

template<bool Local, int Mem, typename MatchFunc, typename StringType>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_dispatch(const StringType& seq1, const StringType& seq2, uint64_t max_mem) const {
    
    if (seq1.size() + seq2.size() > std::numeric_limits<int32_t>::max()) {
        std::stringstream strm;
        strm << "error:[WFAligner] WFA implementation can only align sequences of combined length <= "  << std::numeric_limits<int32_t>::max() << '\n';
        throw std::runtime_error(strm.str());
    }
    
    // make a match function as directed by the calling environment
    MatchFunc match_func(seq1, seq2);
    
    // decide how small of integers we can get away using, then do either the
    // low memory or standard algorithm with the desired configuration
    std::pair<std::vector<CIGAROp>, int32_t> result;
    size_t max_anti_diag = seq1.size() + seq2.size();
    if (max_anti_diag <= (size_t) std::numeric_limits<int8_t>::max()) {
        switch (Mem) {
            case AdaptiveMem:
                result = wavefront_align_core<Local, true, int8_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case StdMem:
                result = wavefront_align_core<Local, false, int8_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case LowMem:
                result = wavefront_align_low_mem_core<Local, int8_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case Recursive:
                result = wavefront_align_recursive_core<Local, int8_t>(seq1, seq2, match_func);
                break;
        }
    }
    else if (max_anti_diag <= (size_t) std::numeric_limits<int16_t>::max()) {
        switch (Mem) {
            case AdaptiveMem:
                result = wavefront_align_core<Local, true, int16_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case StdMem:
                result = wavefront_align_core<Local, false, int16_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case LowMem:
                result = wavefront_align_low_mem_core<Local, int16_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case Recursive:
                result = wavefront_align_recursive_core<Local, int16_t>(seq1, seq2, match_func);
                break;
        }
    }
    else {
        switch (Mem) {
            case AdaptiveMem:
                result = wavefront_align_core<Local, true, int32_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case StdMem:
                result = wavefront_align_core<Local, false, int32_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case LowMem:
                result = wavefront_align_low_mem_core<Local, int32_t>(seq1, seq2, match_func, max_mem);
                break;
                
            case Recursive:
                result = wavefront_align_recursive_core<Local, int32_t>(seq1, seq2, match_func);
                break;
        }
    }
    // TODO: not doing int64_t for now -- we have bigger problems if we're aligning
    // gigabase sequences
    
    // unreduce the score
    result.second *= factor;
    
    // for local alignment, the conversion back will take place in the local core
    // because we need to measure the length of the aligned portion
    if (match && !Local) {
        result.second = convert_score(seq1.size(), seq2.size(), result.second);
    }
    
    return result;
}

template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_core(const StringType& seq1, const StringType& seq2,
                                const MatchFunc& match_func, uint64_t max_mem) const {

    // the return value
    std::pair<std::vector<CIGAROp>, int32_t> res;
    
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    
    uint64_t curr_mem = 0;
    bool hit_max_mem = false;
    
    // init wavefront to the upper left of both sequences' starts
    std::vector<Wavefront<IntType>> wfs;
    wfs.emplace_back(0, 1);
    wfs.front().M[0] = 0;
    
    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wfs.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, 0, wfs.front(),
                       opt, opt_diag, opt_s, max_s);
    }
    if (Adaptive) {
        curr_mem += wfs.front().bytes();
    }
    while (!wavefront_reached(wfs.back(), final_diag, final_anti_diag) &&
           wfs.size() - 1 <= max_s) {
        if (Adaptive && curr_mem > max_mem) {
            hit_max_mem = true;
            break;
        }
        
        wfs.emplace_back(wavefront_next<IntType>(seq1, seq2, wfs));
        
        wavefront_extend(seq1, seq2, wfs.back(), match_func);
        if (Local) {
            find_local_opt(seq1, seq2, wfs.size() - 1, wfs.back(),
                           opt, opt_diag, opt_s, max_s);
        }
        // prune lagging diagonals
        wavefront_prune<Local>(wfs.back());
        
        if (Adaptive) {
            curr_mem += wfs.back().bytes();
        }
    }
    
#ifdef debug_viz
    wavefront_viz(seq1, seq2, wfs, scores);
#endif
    
    if (Adaptive && hit_max_mem) {
        // we hit the memory maximum without completing, fall back to the low memory algorithm
        res = fall_back_to_low_mem<Local, IntType>(seq1, seq2, match_func, max_mem,
                                                   opt, opt_diag, opt_s, max_s, wfs);
    }
    else {
        // we completed the algorithm
        int64_t traceback_s, traceback_d;
        if (Local) {
            traceback_d = opt_diag;
            traceback_s = opt_s;
        }
        else {
            traceback_d = seq1.size() - seq2.size();
            traceback_s = wfs.size() - 1;
        }
        res.second = traceback_s;
        res.first = wavefront_traceback(seq1, seq2, wfs, traceback_s, traceback_d);
    }
    
    return res;
}

template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::fall_back_to_low_mem(const StringType& seq1, const StringType& seq2,
                                const MatchFunc& match_func, uint64_t max_mem,
                                int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                std::vector<WFAligner::Wavefront<IntType>>& wfs) const {
    
    std::deque<Wavefront<IntType>> wf_buffer;
    std::vector<std::pair<int32_t, Wavefront<IntType>>> wf_bank;
    
    // subsampling params
    int64_t epoch_len = 1;
    int64_t epoch_end = 1;
    int64_t sample_rate = 1;
    
    uint64_t curr_mem_bank = 0;
    uint64_t curr_mem_block = 0;
    uint64_t max_mem_block = 0;
    uint64_t curr_mem_buffer = 0;
    
    size_t s = wfs.size() - 1;
    if (wfs.size() < 2 * stripe_width) {
        // we wouldn't have started evicting yet, everything should end up in the buffer
        for (auto& wf : wfs) {
            curr_mem_buffer += wf.bytes();
            wf_buffer.emplace_back(std::move(wf));
        }
    }
    else {
        // grab the stripe(s) that would have been found in the buffer
        size_t buffer_stripe_num = (wfs.size() + 1) / stripe_width - 1;
        for (size_t i = buffer_stripe_num * stripe_width; i < wfs.size(); ++i) {
            wf_buffer.emplace_back(std::move(wfs[i]));
        }
        
        // sub-sample the stripes according the the epoch schedule
        // TODO: kind of a lot of copy-pasta code here...
        for (int64_t stripe_num = 0; stripe_num < buffer_stripe_num; ++stripe_num) {
            
            // is this the first wavefront of a new epoch?
            if (stripe_num == epoch_end) {
                // we're entering a new epoch, adust the sub-sampling parameters
                // accordingly
                epoch_len *= 4;
                epoch_end += epoch_len;
                sample_rate *= 2;
            }
            
            // bit-hacky sub for % that works because sample_rate is a power of 2
            bool keep = !stripe_num || (stripe_num & (sample_rate - 1));
            for (size_t i = 0, k = stripe_num * stripe_width; i < stripe_width; ++i) {
                if (keep) {
                    // we're ejecting a stripe that's being retained,
                    wf_bank.emplace_back(k + i, std::move(wfs[k + i]));
                    curr_mem_bank += wf_bank.back().second.bytes();
                    curr_mem_block = 0;
                }
                else {
                    curr_mem_block += wfs[k + i].bytes();
                    max_mem_block = std::max(max_mem_block, curr_mem_block);
                }
            }
        }
    }
    
    // clear the memory that we're no longer using
    wfs = std::vector<Wavefront<IntType>>();
    
    return wavefront_align_low_mem_core_internal<Local, true, IntType>(seq1, seq2, match_func, max_mem,
                                                                       epoch_len, epoch_end, sample_rate,
                                                                       opt, opt_diag, opt_s, max_s,
                                                                       wf_buffer, s, wf_bank,
                                                                       curr_mem_buffer, curr_mem_bank, curr_mem_block, max_mem_block);
}

template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_low_mem_core(const StringType& seq1, const StringType& seq2,
                                        const MatchFunc& match_func, uint64_t max_mem) const {

    // subsampling params
    int64_t epoch_len = 1;
    int64_t epoch_end = 1;
    int64_t sample_rate = 1;
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    
    // init wavefront to the upper left of both sequences' starts
    std::deque<Wavefront<IntType>> wf_buffer;
    wf_buffer.emplace_back(0, 1);
    wf_buffer.front().M[0] = 0;
    size_t s = 0;
    
    // make sure the first wavefront is extended
    wavefront_extend(seq1, seq2, wf_buffer.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, s, wf_buffer.front(),
                       opt, opt_diag, opt_s, max_s);
    }
    
    // the wavefronts we've decided to keep after they leave the buffer
    std::vector<std::pair<int32_t, Wavefront<IntType>>> wf_bank;
    
    uint64_t curr_mem = wf_buffer.front().bytes();
    
    // TODO: i'm not sure if all of the tracking variables will be able to get optimized out
    // since it will be hard to identify them as constant. compiler might not respect inlining
    // of internal method
    return wavefront_align_low_mem_core_internal<Local, false, IntType>(seq1, seq2, match_func, max_mem,
                                                                        epoch_len, epoch_end, sample_rate,
                                                                        opt, opt_diag, opt_s, max_s,
                                                                        wf_buffer, s, wf_bank, curr_mem, 0, 0, 0);
};

template<bool Local, bool Adaptive, typename IntType, typename StringType, typename MatchFunc>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_low_mem_core_internal(const StringType& seq1, const StringType& seq2,
                                                 const MatchFunc& match_func, uint64_t max_mem,
                                                 int64_t epoch_len, int64_t epoch_end, int64_t sample_rate,
                                                 int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                                 std::deque<WFAligner::Wavefront<IntType>>& wf_buffer, size_t s,
                                                 std::vector<std::pair<int32_t, WFAligner::Wavefront<IntType>>>& wf_bank,
                                                 uint64_t curr_mem_buffer, uint64_t curr_mem_bank,
                                                 uint64_t curr_mem_block, uint64_t max_mem_block) const {

    // the return value
    std::pair<std::vector<CIGAROp>, int32_t> res;
    
    // use these to know when we've finished the alignment
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();
    
    // FIXME: i also need an estimate of the recomputed block...
    bool hit_max_mem = false;
    
    // do wavefront iterations until hit max score or alignment finishes
    while (!wavefront_reached(wf_buffer.back(), final_diag, final_anti_diag) &&
           s <= max_s) {
        
        // increase score and get the next wavefront
        ++s;
        wf_buffer.emplace_back(wavefront_next<IntType>(seq1, seq2, wf_buffer));
        if (Adaptive) {
            curr_mem_buffer += wf_buffer.back().bytes();
        }
        
        // note: the method of keeping up to 2 stripes in memory at a time is a bit strange
        // here, but it greatly simplifies the fall back to the recursive alignment
        if (wf_buffer.size() == 2 * stripe_width) {
            // we need to evict a stripe from the buffer to continue
            
            // the stripe of the wavefront that is falling out of the buffer
            int32_t stripe_num = s / stripe_width - 1; // formula works because s must be k*w-1
            
            // is this the first wavefront of a new epoch?
            if (stripe_num == epoch_end) {
                // we're entering a new epoch, adust the sub-sampling parameters
                // accordingly
                epoch_len *= 4;
                epoch_end += epoch_len;
                sample_rate *= 2;
            }
            
            // bit-hacky sub for % that works because sample_rate is a power of 2
            bool keep = !stripe_num || (stripe_num & (sample_rate - 1));
            for (size_t i = 0; i < stripe_width; ++i) {
                if (keep) {
                    // we're ejecting a stripe that's being retained,
                    wf_bank.emplace_back(stripe_num * stripe_width + i, std::move(wf_buffer.front()));
                    if (Adaptive) {
                        curr_mem_block = 0;
                        curr_mem_bank += wf_buffer.front().bytes();
                    }
                }
                else if (Adaptive) {
                    // record that we're freeing memory
                    auto b = wf_buffer.front().bytes();
                    curr_mem_buffer -= b;
                    curr_mem_block += b;
                    max_mem_block = std::max(curr_mem_block, max_mem_block);
                }
                wf_buffer.pop_front();
            }
        }
        
        // follow matches
        wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
        
        if (Local) {
            find_local_opt(seq1, seq2, s, wf_buffer.back(),
                           opt, opt_diag, opt_s, max_s);
        }
        // prune lagging diagonals
        wavefront_prune<Local>(wf_buffer.back());
        
        if (Adaptive) {
            curr_mem_bank += wf_buffer.back().bytes();
            if (curr_mem_bank + max_mem_block + curr_mem_buffer >= max_mem) {
                hit_max_mem = true;
                break;
            }
        }
    }
    
    if (Adaptive && hit_max_mem) {
        // we hit the maximum memory without completing, fall back to the recursive algorithm
        res = fall_back_to_recursive<Local, IntType>(seq1, seq2, match_func,
                                                     opt, opt_diag, opt_s, max_s,
                                                     wf_buffer, s, wf_bank);
    }
    else {
        // we completed successfully
        
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
        
        res.first = wavefront_traceback_low_mem<Local>(seq1, seq2, wf_bank, traceback_s, traceback_d, match_func);
        res.second = s;
    }
    
    return res;
}

template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::fall_back_to_recursive(const StringType& seq1, const StringType& seq2,
                                  const MatchFunc& match_func,
                                  int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                  std::deque<WFAligner::Wavefront<IntType>>& wf_buffer, size_t s,
                                  std::vector<std::pair<int32_t, WFAligner::Wavefront<IntType>>>& wf_bank) const {
    
    // convert the low memory data structures into the recursive data structures
    
    std::vector<Wavefront<IntType>> opt_stripe;
    std::vector<Wavefront<IntType>> first_stripe;
    
    if (!wf_bank.empty()) {
        // we ejected a stripe from the buffer, which we guarantee must include
        // the first stripe, so now we remember it
        first_stripe.reserve(stripe_width);
        for (size_t i = 0; i < stripe_width; ++i) {
            first_stripe.emplace_back(std::move(wf_bank[i].second));
        }
        
        if (Local) {
            
            int64_t opt_stripe_num = (opt_s + stripe_width - 1) / stripe_width;
            if (opt_stripe_num != 0 && opt_stripe_num < (s - wf_buffer.size()) / stripe_width + 1) {
                // the opt isn't in the first stripe or the buffer, so we need to either find or recompute it
                
                size_t i = 0;
                while (wf_bank[i].first / stripe_width < opt_stripe_num) {
                    i += stripe_width;
                }
                
                if (wf_bank[i].first / stripe_width != opt_stripe_num) {
                    // we didn't sub-sample the stripe with the opt, so we have recompute it
                    
                    // eliminate the portion after the opt, which we don't need
                    wf_bank.resize(i);
                    // and add another stripe to the bank
                    WFBankAdapter<IntType> adapter(wf_bank);
                    for (size_t j = 0; j < stripe_width; ++j) {
                        wf_bank.emplace_back(wf_bank.back().first + 1,
                                             wavefront_next<IntType>(seq1, seq2, adapter));
                        wavefront_extend(seq1, seq2, wf_bank.back().second, match_func);
                        wavefront_prune<Local>(wf_bank.back().second);
                    }
                }
                
                // retrieve the opt's stripe from the bank
                opt_stripe.reserve(stripe_width);
                for (size_t n = i + stripe_width; i < n; ++i) {
                    opt_stripe.emplace_back(std::move(wf_bank[i].second));
                }
            }
        }
    }
    
    // discard the wavefronts we're no longer retaining
    wf_bank = std::vector<std::pair<int32_t, Wavefront<IntType>>>();
    
    return wavefront_align_recursive_core_internal<Local, IntType>(seq1, seq2, match_func, opt, opt_diag, opt_s, max_s,
                                                                   opt_stripe, first_stripe, wf_buffer, s);
}

template<bool Local, typename IntType, typename StringType, typename MatchFunc>
std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_recursive_core(const StringType& seq1, const StringType& seq2,
                                          const MatchFunc& match_func) const {
    
    // trackers used for local alignment
    int64_t opt = 0;
    int64_t opt_diag = 0;
    int64_t opt_s = 0;
    size_t max_s = std::numeric_limits<size_t>::max();
    std::vector<Wavefront<IntType>> opt_stripe;
    
    // init wavefront to the upper left of both sequences' starts
    std::deque<Wavefront<IntType>> wf_buffer;
    wf_buffer.emplace_back(0, 1);
    wf_buffer.front().M[0] = 0;
    size_t s = 0;
    
    wavefront_extend(seq1, seq2, wf_buffer.front(), match_func);
    if (Local) {
        find_local_opt(seq1, seq2, s, wf_buffer.front(),
                       opt, opt_diag, opt_s, max_s);
    }
    
    // we need to remember the first stripe for the recursive step
    std::vector<Wavefront<IntType>> first_stripe;
    
    return wavefront_align_recursive_core_internal<Local, IntType>(seq1, seq2, match_func, opt, opt_diag, opt_s, max_s,
                                                                   opt_stripe, first_stripe, wf_buffer, s);
}

template<bool Local, typename IntType, typename StringType, typename MatchFunc>
inline std::pair<std::vector<CIGAROp>, int32_t>
WFAligner::wavefront_align_recursive_core_internal(const StringType& seq1, const StringType& seq2,
                                                   const MatchFunc& match_func,
                                                   int64_t opt, int64_t opt_diag, int64_t opt_s, size_t max_s,
                                                   std::vector<WFAligner::Wavefront<IntType>>& opt_stripe,
                                                   std::vector<WFAligner::Wavefront<IntType>>& first_stripe,
                                                   std::deque<WFAligner::Wavefront<IntType>>& wf_buffer, size_t s) const {
    
    // use these to know when we've finished the alignment
    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();
    
    // do wavefront iterations until hit max score or alignment finishes
    while (!wavefront_reached(wf_buffer.back(), final_diag, final_anti_diag) &&
           s <= max_s) {
        // increase score and get the next wavefront
        ++s;
        wf_buffer.emplace_back(wavefront_next<IntType>(seq1, seq2, wf_buffer));
        
        if (wf_buffer.size() == 2 * stripe_width) {
            // we need to evict a stripe wavefront from the buffer to continue
            for (int64_t i = 0; i < stripe_width; ++i) {
                if (s + 1 == 2 * stripe_width) {
                    // we need to remember the first stripe when we evict it
                    first_stripe.emplace_back(std::move(wf_buffer.front()));
                }
                else if (Local && opt_s > s - 2 * stripe_width && opt_s <= s - stripe_width) {
                    // we're ejecting the stripe that contains the local opt, so
                    // we remember it for later
                    if (i == 0) {
                        // we might have an existing stripe in the opt stripe
                        opt_stripe.clear();
                    }
                    opt_stripe.emplace_back(std::move(wf_buffer.front()));
                }
                wf_buffer.pop_front();
            }
            
        }
        
        // follow matches
        wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
        
        if (Local) {
            find_local_opt(seq1, seq2, s, wf_buffer.back(),
                           opt, opt_diag, opt_s, max_s);
        }
        // prune lagging diagonals
        wavefront_prune<Local>(wf_buffer.back());
    }
    
    int64_t traceback_d;
    int64_t traceback_s;
    int32_t score;
    if (Local) {
        traceback_d = opt_diag;
        traceback_s = opt_s;
        score = opt_s;
        
        // TODO: i don't like repeating this formula 3x in different places...
        int64_t opt_stripe_num = (opt_s + stripe_width - 1) / stripe_width;
        if (opt_stripe_num != 0 && opt_stripe_num != (s - wf_buffer.size()) / stripe_width + 1) {
            // the local optimum did not occur in the first or last stripe, so we replace the contents
            // of the last stripe with the opt's stripe, which we should have retained
            wf_buffer.clear();
            for (auto& wf : opt_stripe) {
                wf_buffer.emplace_back(std::move(wf));
            }
        }
    }
    else {
        traceback_d = final_diag;
        traceback_s = s;
        score = s;
    }
    
    std::vector<CIGAROp> cigar;
    if (s < 2 * stripe_width) {
        // we finished the entire alignment without evicting, so we just do normal traceback
        // TODO: i should be able to handle this with a template, but it's having trouble with
        // type inference...
        std::vector<Wavefront<IntType>> traceback_wfs;
        traceback_wfs.reserve(wf_buffer.size());
        for (auto& wf : wf_buffer) {
            traceback_wfs.emplace_back(std::move(wf));
        }
        cigar = wavefront_traceback(seq1, seq2, traceback_wfs,
                                    traceback_s, traceback_d);
    }
    else if (Local && traceback_s < stripe_width) {
        // we encountered the local optimum within the first stripe, so we can skip complicated
        // traceback logic
        cigar = wavefront_traceback(seq1, seq2, first_stripe,
                                    traceback_s, traceback_d);
    }
    else {
        // we enter the recursive portion of the algorithm
        int64_t stripe_num = (traceback_s - wf_buffer.size()) / stripe_width + 1;
        WFMatrix_t traceback_mat = MAT_M;
        int64_t lead_matches = 0;
        wavefront_align_recursive_internal<Local, IntType>(seq1, seq2,
                                                           0, first_stripe,
                                                           stripe_num, wf_buffer,
                                                           traceback_s, traceback_d,
                                                           traceback_mat, lead_matches,
                                                           cigar, match_func);
        assert(traceback_s == 0 && traceback_d == 0 && traceback_mat == MAT_M);
        // the CIGAR is constructed in reverse, unreverse it
        std::reverse(cigar.begin(), cigar.end());
        // also merge operations from adjacent blocks
        coalesce_cigar(cigar);
    }
    
    
    return std::make_pair(cigar, int32_t(score));
}


template<bool Local, typename IntType, typename StringType, typename MatchFunc, typename WFVector1, typename WFVector2>
void WFAligner::wavefront_align_recursive_internal(const StringType& seq1, const StringType& seq2,
                                                   int64_t lower_stripe_num, WFVector1& lower_stripe,
                                                   int64_t upper_stripe_num, WFVector2& upper_stripe,
                                                   int64_t& traceback_s, int64_t& traceback_d,
                                                   WFMatrix_t& traceback_mat, int64_t& lead_matches,
                                                   std::vector<CIGAROp>& cigar,
                                                   const MatchFunc& match_func) const {
    
    if (lower_stripe_num + 1 == upper_stripe_num) {
        // do a step of traceback
        
        // collect all of the wavefronts in one vector
        for (auto& wf : upper_stripe) {
            lower_stripe.emplace_back(std::move(wf));
        }
        
        
        // convert the current s into an index within this block
        int64_t relative_s = traceback_s - lower_stripe_num * stripe_width;
        
        // traceback the stripe
        wavefront_traceback_internal(seq1, seq2, lower_stripe,
                                     traceback_d, lead_matches, traceback_mat,
                                     relative_s, cigar);
        
        // convert the final within-block s back into the real s
        traceback_s = lower_stripe_num * stripe_width + relative_s;
        
        // we no longer need the upper stripe for anything, we can discard it
        while (lower_stripe.size() > stripe_width) {
            lower_stripe.pop_back();
        }
    }
    else {
        // borrow the parent call's wavefronts to do our own DP with
        std::deque<Wavefront<IntType>> wf_buffer;
        //wf_buffer.reserve(2 * stripe_width);
        for (auto& wf : lower_stripe) {
            wf_buffer.emplace_back(std::move(wf));
        }
        
        // do DP out to the middle of the row interval
        int64_t target_stripe_num = (lower_stripe_num + upper_stripe_num) / 2;
        for (int64_t i = 0, n = stripe_width * (target_stripe_num - lower_stripe_num); i < n; ++i) {
            // one step of DP
            wf_buffer.emplace_back(wavefront_next<IntType>(seq1, seq2, wf_buffer));
            wavefront_extend(seq1, seq2, wf_buffer.back(), match_func);
            // prune lagging diagonals
            // TODO: need to update this to Local when it's implemented
            wavefront_prune<Local>(wf_buffer.back());
            
            if (wf_buffer.size() == 2 * stripe_width) {
                // we've completed a stripe, we can now evict the previous one
                for (int64_t j = 0; j < stripe_width; ++j) {
                    if (i < stripe_width) {
                        // we're evicting the first stripe, so we have to give it back
                        // to the parent call
                        lower_stripe[j] = std::move(wf_buffer.front());
                        
                    }
                    wf_buffer.pop_front();
                }
            }
        }
        
        // recurse into upper half
        wavefront_align_recursive_internal<Local, IntType>(seq1, seq2,
                                                           target_stripe_num, wf_buffer,
                                                           upper_stripe_num, upper_stripe,
                                                           traceback_s, traceback_d,
                                                           traceback_mat, lead_matches,
                                                           cigar, match_func);
        // recurse into lower half
        wavefront_align_recursive_internal<Local, IntType>(seq1, seq2,
                                                           lower_stripe_num, lower_stripe,
                                                           target_stripe_num, wf_buffer,
                                                           traceback_s, traceback_d,
                                                           traceback_mat, lead_matches,
                                                           cigar, match_func);
    }
}

template<typename PrefSemilocalWFA, typename AnchorGlobalWFA, typename SuffSemilocalWFA>
std::tuple<std::vector<CIGAROp>, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>>
WFAligner::wavefront_align_local_core(const char* seq1, size_t len1,
                                      const char* seq2, size_t len2,
                                      uint64_t max_mem,
                                      size_t anchor_begin_1, size_t anchor_end_1,
                                      size_t anchor_begin_2, size_t anchor_end_2,
                                      bool anchor_is_match) const {
    
    if (match == 0) {
        std::stringstream strm;
        strm << "error:[WFAligner] aligner must be initialized with Smith-Waterman-Gotoh-style scoring parameters to perform local alignment" << std::endl;
        throw std::runtime_error(strm.str());
    }
    
    // tail align the prefix
    RevStringView pref1(seq1, 0, anchor_begin_1);
    RevStringView pref2(seq2, 0, anchor_begin_2);
    
    auto pref_result = PrefSemilocalWFA(this, max_mem)(pref1, pref2);
    
    std::vector<CIGAROp> cigar = std::move(pref_result.first);
    std::reverse(cigar.begin(), cigar.end());
    auto aligned_pref_len = cigar_base_length(cigar);
    int32_t pref_score = convert_score(aligned_pref_len.first, aligned_pref_len.second,
                                       pref_result.second);
    
    int32_t anchor_score = 0;
    if (!anchor_is_match) {
        
        // align the anchor
        StringView anchor1(seq1, anchor_begin_1, anchor_end_1 - anchor_begin_1);
        StringView anchor2(seq2, anchor_begin_2, anchor_end_2 - anchor_begin_2);
        
        auto anchor_res = AnchorGlobalWFA(this, max_mem)(anchor1, anchor2);
        
        // get the score and add to the CIGAR
        anchor_score = convert_score(anchor_end_1 - anchor_begin_1, anchor_end_2 - anchor_begin_2,
                                     anchor_res.second);
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
        anchor_score = match * len;
        cigar.emplace_back('M', len);
    }
    
    // tail align the suffix
    StringView suff1(seq1, anchor_end_1, len1 - anchor_end_1);
    StringView suff2(seq2, anchor_end_2, len2 - anchor_end_2);
    
    auto suff_result = SuffSemilocalWFA(this, max_mem)(suff1, suff2);
    
    auto aligned_suff_len = cigar_base_length(suff_result.first);
    int32_t suff_score = convert_score(aligned_suff_len.first, aligned_suff_len.second,
                                       suff_result.second);
    
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

template<typename StringType, typename MatchFunc, typename IntType>
inline void WFAligner::wavefront_extend(const StringType& seq1, const StringType& seq2,
                                        Wavefront<IntType>& wf, const MatchFunc& match_func) const {
    
    for (int64_t k = 0; k < wf.size(); ++k) {
        int64_t diag = wf.diag_begin + k;
        int64_t anti_diag = wf.M[k];
        int64_t i = (diag + anti_diag) / 2;
        int64_t j = (anti_diag - diag) / 2;
        if (i >= 0 && j >= 0) {
            wf.M[k] += 2 * match_func(i, j);
        }
    }
}

template<typename IntType, typename WFVector, typename StringType>
inline WFAligner::Wavefront<IntType> WFAligner::wavefront_next(const StringType& seq1, const StringType& seq2,
                                                               const WFVector& wfs) const {
    
    // find the range of incoming wavefronts
    int32_t hi = std::numeric_limits<int32_t>::min();
    int32_t lo = std::numeric_limits<int32_t>::max();
    // mismatch stays in the same diagonals
    if (wfs.size() >= mismatch) {
        const auto& wf_prev = wfs[wfs.size() - mismatch];
        if (!wf_prev.empty()) {
            lo = std::min<int32_t>(lo, wf_prev.diag_begin);
            hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.size());
        }
    }
    // TODO: but sometimes these bounds are slightly too large because of the gap extend
    // but it should be only ~s * o extra work
    // insertion/deletion expand the diagonals by 1
    for (auto s_back : {gap_extend, gap_open + gap_extend}) {
        if (wfs.size() >= s_back) {
            const auto& wf_prev = wfs[wfs.size() - s_back];
            if (!wf_prev.empty()) {
                lo = std::min<int32_t>(lo, wf_prev.diag_begin - 1);
                hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.size() + 1);
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
    
    Wavefront<IntType> wf(ctor_lo, ctor_hi);
    
    if (lo < hi) {
        
        // open gaps
        if (wfs.size() >= gap_extend + gap_open) {
            const auto& wf_prev = wfs[wfs.size() - gap_extend - gap_open];
            int64_t prev_diag_end = wf_prev.diag_begin + wf_prev.size();
            for (auto k = lo; k < hi; ++k) {
                if (k - 1 >= wf_prev.diag_begin && k - 1 < prev_diag_end) {
                    int64_t a = wf_prev.M[k - 1 - wf_prev.diag_begin] + 1;
                    if ((a + k) / 2 <= (int64_t) seq1.size() && (a - k) / 2 <= (int64_t) seq2.size()) {
                        wf.I[k - lo] = a;
                    }
                }
                if (k + 1 >= wf_prev.diag_begin && k + 1 < prev_diag_end) {
                    int64_t a = wf_prev.M[k + 1 - wf_prev.diag_begin] + 1;
                    if ((a + k) / 2 <= (int64_t) seq1.size() && (a - k) / 2 <= (int64_t) seq2.size()) {
                        wf.D[k - lo] = a;
                    }
                }
            }
        }
        
        // extend gaps
        if (wfs.size() >= gap_extend) {
            const auto& wf_prev = wfs[wfs.size() - gap_extend];
            int64_t prev_diag_end = wf_prev.diag_begin + wf_prev.size();
            for (auto k = lo; k < hi; ++k) {
                if (k - 1 >= wf_prev.diag_begin && k - 1 < prev_diag_end) {
                    int64_t a = wf_prev.I[k - 1 - wf_prev.diag_begin] + 1;
                    if ((a + k) / 2 <= (int64_t) seq1.size() && (a - k) / 2 <= (int64_t) seq2.size()) {
                        wf.I[k - lo] = std::max<int32_t>(wf.I[k - lo], a);
                    }
                }
                if (k + 1 >= wf_prev.diag_begin && k + 1 < prev_diag_end) {
                    int64_t a = wf_prev.D[k + 1 - wf_prev.diag_begin] + 1;
                    if ((a + k) / 2 <= (int64_t) seq1.size() && (a - k) / 2 <= (int64_t) seq2.size()) {
                        wf.D[k - lo] = std::max<int32_t>(wf.D[k - lo], a);
                    }
                }
            }
        }
        
        // mismatches
        if (wfs.size() >= mismatch) {
            const auto& wf_prev = wfs[wfs.size() - mismatch];
            int64_t prev_diag_end = wf_prev.diag_begin + wf_prev.size();
            for (auto k = lo; k < hi; ++k) {
                if (k >= wf_prev.diag_begin && k < prev_diag_end) {
                    int64_t a = wf_prev.M[k - wf_prev.diag_begin] + 2;
                    if ((a + k) / 2 <= (int64_t) seq1.size() && (a - k) / 2 <= (int64_t) seq2.size()) {
                        wf.M[k - lo] = a;
                    }
                }
            }
        }
        
        // calculate the opt
        for (size_t k = 0; k < wf.size(); ++k) {
            wf.M[k] = std::max(wf.M[k], std::max(wf.I[k], wf.D[k]));
        }
    }
    
    return wf;
}

template<bool Local, typename IntType>
inline void WFAligner::wavefront_prune(WFAligner::Wavefront<IntType>& wf) const {
    // TODO: for now, just skipping pruning for local alignment, but actually should come up
    // with a better pruning mechanism for it...
    if (!wf.empty() && lagging_diagonal_prune >= 0 && !Local) {
        IntType max_anti_diag = std::numeric_limits<IntType>::min();
        for (int32_t i = 0; i < wf.size(); ++i) {
            max_anti_diag = std::max(max_anti_diag, wf.M[i]);
        }
        
        while (wf.M[0] < max_anti_diag - lagging_diagonal_prune) {
            wf.pop_front();
        }
        
        while (wf.M[wf.len - 1] < max_anti_diag - lagging_diagonal_prune) {
            wf.pop_back();
        }
    }
}

template<typename IntType>
inline bool WFAligner::wavefront_reached(const WFAligner::Wavefront<IntType>& wf,
                                         int32_t diag, int32_t anti_diag) const {
    if (diag >= wf.diag_begin && diag < wf.diag_begin + wf.size()) {
        return (wf.M[diag - wf.diag_begin] == anti_diag);
    }
    return false;
}

template <typename StringType, typename IntType>
inline std::vector<CIGAROp> WFAligner::wavefront_traceback(const StringType& seq1, const StringType& seq2,
                                                           const std::vector<WFAligner::Wavefront<IntType>>& wfs,
                                                           int64_t s, int64_t d) const {
    
    // always begin traceback in the match matrix
    WFMatrix_t mat = MAT_M;
    int64_t lead_matches = 0;
    
    std::vector<CIGAROp> cigar;
    wavefront_traceback_internal(seq1, seq2, wfs, d, lead_matches, mat, s, cigar);
    assert(s == 0 && d == 0 && mat == MAT_M);
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

template <bool Local, typename StringType, typename MatchFunc, typename IntType>
std::vector<CIGAROp> WFAligner::wavefront_traceback_low_mem(const StringType& seq1, const StringType& seq2,
                                                            std::vector<std::pair<int32_t, Wavefront<IntType>>>& wf_bank,
                                                            int64_t s, int64_t d, const MatchFunc& match_func) const {
    
    // the return value
    std::vector<CIGAROp> cigar;
    
    if (Local) {
        // if we overshot the opt traceback location, clear the way for it to be
        // in the last position of the bank
        while (wf_bank.back().first > s) {
            wf_bank.pop_back();
        }
        // and recompute out to the opt if necessary
        WFBankAdapter<IntType> adapter(wf_bank);
        while (wf_bank.back().first < s) {
            wf_bank.emplace_back(wf_bank.back().first + 1,
                                 wavefront_next<IntType>(seq1, seq2, adapter));
            wavefront_extend(seq1, seq2, wf_bank.back().second, match_func);
            wavefront_prune<Local>(wf_bank.back().second);
        }
        // the traceback starting location should now be in the final wavefront
        // in the bank, regardless of global/local
    }
    
    // a moving window over the DP that we will incrementally do traceback through
    std::deque<Wavefront<IntType>> traceback_block;
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
    int64_t relative_s = traceback_block.size() - 1;
    
    // traceback as far as possible in this block
    wavefront_traceback_internal(seq1, seq2, traceback_block,
                                 d, lead_matches, mat, relative_s, cigar);
    // remove anything we don't need anymore
    traceback_block.resize(relative_s + 1);
    
    while (i > 0) {
        // there's an earlier discontiguous stripe we need to connect to
        
        // figure out where it starts
        int64_t j = i - 1;
        while (j - 1 >= 0 && wf_bank[j - 1].first == wf_bank[j].first - 1) {
            --j;
        }
        
        // initialize the full block with stripe
        std::vector<Wavefront<IntType>> fillin_wfs;
        fillin_wfs.reserve(wf_bank[i].first - wf_bank[j].first);
        for (int64_t k = j; k < i; ++k) {
            fillin_wfs.emplace_back(std::move(wf_bank[k].second));
        }
        
        // redo the DP between the stripes
        for (int32_t t = wf_bank[i - 1].first + 1, end = wf_bank[i].first; t < end; ++t) {
            fillin_wfs.emplace_back(wavefront_next<IntType>(seq1, seq2, fillin_wfs));
            wavefront_extend(seq1, seq2, fillin_wfs.back(), match_func);
            wavefront_prune<Local>(fillin_wfs.back());
        }
        
        // and move the redone wavefronts into the block buffer for traceback
        for (int64_t k = fillin_wfs.size() - 1; k >= 0; --k) {
            traceback_block.emplace_front(std::move(fillin_wfs[k]));
        }
        i = j;
        
        // traceback as far as possible in this block
        relative_s = traceback_block.size() - 1;
        wavefront_traceback_internal(seq1, seq2, traceback_block,
                                     d, lead_matches, mat, relative_s, cigar);
        // remove anything we don't need anymore
        traceback_block.resize(relative_s + 1);
    }
    
    // we sometimes need to coalesce CIGAR operations that occurred across
    // the boundaries of blocks
    coalesce_cigar(cigar);
    // the CIGAR is built backwards
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

template<typename WFVector, typename StringType>
void WFAligner::wavefront_traceback_internal(const StringType& seq1, const StringType& seq2,
                                             const WFVector& wfs, int64_t& d, int64_t& lead_matches,
                                             WFMatrix_t& mat, int64_t& s, std::vector<CIGAROp>& cigar) const {
    uint32_t op_len = 0;
    while (true) {
        // we're in the match/mismatch matrix
        const auto& wf = wfs[s];
        if (mat == MAT_M) {
            int64_t a = wf.M[d - wf.diag_begin];
            // TODO: i don't love this solution for the partial traceback problem...
            a -= 2 * lead_matches;
            lead_matches = 0;
            // extend through any matches
            int64_t i = (d + a) / 2 - 1;
            int64_t j = (a - d) / 2 - 1;
            while (i >= 0 && j >= 0 && seq1[i] == seq2[j]
                   && a != wf.I[d - wf.diag_begin] && a != wf.D[d - wf.diag_begin]) {
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
            if (a == wf.I[d - wf.diag_begin]) {
                // this is where an insertion closed
                if (op_len != 0) {
                    cigar.emplace_back('M', op_len);
                }
                mat = MAT_I;
                op_len = 0;
                lead_matches = 0;
                continue;
            }
            else if (a == wf.D[d - wf.diag_begin]) {
                // this is where a deletion closed
                if (op_len != 0) {
                    cigar.emplace_back('M', op_len);
                }
                mat = MAT_D;
                op_len = 0;
                lead_matches = 0;
                continue;
            }
            else if (s >= mismatch) {
                // this must a mismatch
                op_len += 1;
                s -= mismatch;
                lead_matches = 0;
                continue;
            }
            // the traceback goes outside of the DP structure (should only happen
            // in the low-memory code path)
            break;
        }
        else if (mat == MAT_I) {
            // we're in the insertion matrix
            if (s >= gap_extend) {
                const auto& wf_prev = wfs[s - gap_extend];
                if (d - 1 >= wf_prev.diag_begin && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()) {
                    if (wf.I[d - wf.diag_begin] == wf_prev.I[d - 1 - wf_prev.diag_begin] + 1) {
                        // an insert extended here
                        s -= gap_extend;
                        op_len += 1;
                        d -= 1;
                        continue;
                    }
                }
            }
            if (s >= gap_extend + gap_open) {
                const auto& wf_prev = wfs[s - gap_extend - gap_open];
                if (d - 1 >= wf_prev.diag_begin  && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()) {
                    if (wf.I[d - wf.diag_begin] == wf_prev.M[d - 1 - wf_prev.diag_begin] + 1) {
                        // an insert opened here
                        s -= gap_extend + gap_open;
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
            if (s >= gap_extend) {
                const auto& wf_prev = wfs[s - gap_extend];
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()) {
                    if (wf.D[d - wf.diag_begin] == wf_prev.D[d + 1 - wf_prev.diag_begin] + 1) {
                        // a deletion extended here
                        s -= gap_extend;
                        op_len += 1;
                        d += 1;
                        continue;
                    }
                }
            }
            if (s >= gap_extend + gap_open) {
                const auto& wf_prev = wfs[s - gap_extend - gap_open];
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.size()){
                    if (wf.D[d - wf.diag_begin] == wf_prev.M[d + 1 - wf_prev.diag_begin] + 1) {
                        // a deletion opened here
                        s -= gap_extend + gap_open;
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


uint32_t WFAligner::gcd(uint32_t a, uint32_t b) const {
    // euclid's algorithm
    if (a < b) {
        std::swap(a, b);
    }
    auto r = a % b;
    if (r == 0) {
        return b;
    }
    else {
        return gcd(b, r);
    }
}

template<typename StringType, typename IntType>
inline void WFAligner::find_local_opt(const StringType& seq1, const StringType& seq2,
                                       int64_t s, const Wavefront<IntType>& wf,
                                       int64_t& opt, int64_t& opt_diag, int64_t& opt_s, size_t& max_s) const {
    
    for (int64_t k = 0; k < wf.size(); ++k) {
        // the score of the corresponding local alignment
        auto a = wf.M[k];
        int32_t local_s = match * a - s;
        if (local_s > opt && a >= std::abs(wf.diag_begin + k)) {
            // we found a new best that is inside the matrix
            opt = local_s;
            opt_diag = wf.diag_begin + k;
            opt_s = s;
            // update our bound on how far out we need to look
            max_s = s + match * (2 * std::min(seq1.size(), seq2.size()) - a);
        }
    }
}

inline int32_t WFAligner::convert_score(size_t len1, size_t len2, int32_t score) const {
    return (match * (len1 + len2) - score) / 2;
}

inline std::pair<size_t, size_t> WFAligner::cigar_base_length(const std::vector<CIGAROp>& cigar) const {
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
                break;
            case 'D':
                len2 += c.len;
                break;
                
            default:
                break;
        }
    }
    return std::make_pair(len1, len2);
}

// merge all adjacent, equivalent operations
inline void WFAligner::coalesce_cigar(std::vector<CIGAROp>& cigar) const {
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

template<typename IntType>
inline void WFAligner::print_wf(const WFAligner::Wavefront<IntType>& wf, int diag_begin, int diag_end) const {
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.size()) {
            int a = wf.M[d - wf.diag_begin];
            if (a > -100) {
                std::cerr << (int) a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.size()) {
            int a = wf.I[d - wf.diag_begin];
            if (a > -100) {
                std::cerr << (int) a;
            }
            else {
                std::cerr << ".";
            }
        }
        std::cerr << "\t";
    }
    std::cerr << std::endl;
    for (int d = diag_begin; d < diag_end; ++d) {
        if (d >= wf.diag_begin && d < wf.diag_begin + (int) wf.size()) {
            int a = wf.D[d - wf.diag_begin];
            if (a > -100) {
                std::cerr << (int) a;
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
void WFAligner::print_wfs(const WFVector& wfs) const {
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int s = 0; s < wfs.size(); ++s) {
        min_begin = std::min<int>(wfs[s].diag_begin, min_begin);
        max_end = std::max<int>(wfs[s].diag_begin + wfs[s].size(), max_end);
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

template<typename IntType>
inline void WFAligner::print_wfs_low_mem(const std::vector<std::pair<int32_t, WFAligner::Wavefront<IntType>>>& wf_bank,
                                         const std::deque<WFAligner::Wavefront<IntType>>& wf_buffer, int s) const {
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int i = 0; i < wf_bank.size(); ++i) {
        min_begin = std::min<int>(wf_bank[i].second.diag_begin, min_begin);
        max_end = std::max<int>(wf_bank[i].second.diag_begin + wf_bank[i].second.size(), max_end);
    }
    for (int i = 0; i < wf_buffer.size(); ++i) {
        min_begin = std::min<int>(wf_buffer[i].diag_begin, min_begin);
        max_end = std::max<int>(wf_buffer[i].diag_begin + wf_buffer[i].size(), max_end);
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

template<typename IntType>
inline void WFAligner::print_wfs_tb(const std::deque<WFAligner::Wavefront<IntType>>& traceback_block,
                                    const std::vector<std::pair<int32_t, WFAligner::Wavefront<IntType>>>& wf_bank,
                                    int i) const {

    std::cerr << "### print traceback block " << wf_bank[i].first << ":" << wf_bank[i].first + traceback_block.size() << " ###" << std::endl;
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int k = 0; k < traceback_block.size(); ++k) {
        min_begin = std::min<int>(min_begin, traceback_block[k].diag_begin);
        max_end = std::max<int>(max_end, traceback_block[k].diag_begin + traceback_block[k].size());
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

template <typename StringType, typename WFVector>
void WFAligner::wavefront_viz(const StringType& seq1, const StringType& seq2,
                              const WFVector& wfs) const {

    std::vector<std::vector<int>> matrix(seq1.size() + 1, std::vector<int>(seq2.size() + 1, -1));

    std::vector<std::vector<bool>> opt(seq1.size() + 1, std::vector<bool>(seq2.size() + 1, false));


    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size();

    if (wavefront_reached(wfs[wfs.size() - 1], final_diag, final_anti_diag)) {
        auto cigar = wavefront_traceback(seq1, seq2, wfs, wfs.size() - 1, final_diag);
        int i = 0, j = 0;
        opt[i][j] = true;
        for (auto& op : cigar) {
            for (int k = 0; k < op.len; ++k) {
                if (op.op == 'M' || op.op == 'I') {
                    ++i;
                }
                if (op.op == 'M' || op.op == 'D') {
                    ++j;
                }
                opt[i][j] = true;
            }
        }
    }

    for (int s = 0; s < wfs.size(); ++s) {
        const auto& wf = wfs[s];
        for (int k = 0; k < wf.size(); ++k) {
            int d = wf.diag_begin + k;
            int a = wf.M[k];
            if (a < -2) {
                continue;
            }
            int i = (a + d) / 2 - 1;
            int j = (a - d) / 2 - 1;
            while (i >= 0 && j >= 0 && i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]
                   && a != wf.I[k] && a != wf.D[k]) {
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
            if (opt[i][j]) {
                std::cerr << '*';
            }
        }
        std::cerr << std::endl;
    }
}

} // end namespace wfalm


#endif /* wfa_lm.hpp */
