// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_ALIGN_PARALLEL_DP_THREAD_LOCAL_STORAGE_H_
#define SEQAN_INCLUDE_ALIGN_PARALLEL_DP_THREAD_LOCAL_STORAGE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TConfig>
class WavefrontAlignmentThreadLocalStorage
{
public:
    //-------------------------------------------------------------------------
    // Member Types.

    using TIntermediate = typename TConfig::TIntermediate;
    using TCache        = typename TConfig::TCache;

    using TAlignmentLocal = std::tuple<TIntermediate, TCache>;

    //-------------------------------------------------------------------------
    // Private Members.

    std::vector<TAlignmentLocal>   _mMultiAlignmentThreadLocal;

    //-------------------------------------------------------------------------
    // Constructor.

    explicit WavefrontAlignmentThreadLocalStorage(size_t const numAlignments) :
        _mMultiAlignmentThreadLocal(numAlignments)
    {}

    // Delegating default constructor.
    WavefrontAlignmentThreadLocalStorage() : WavefrontAlignmentThreadLocalStorage(1)
    {}

    WavefrontAlignmentThreadLocalStorage(WavefrontAlignmentThreadLocalStorage const &) = default;
    WavefrontAlignmentThreadLocalStorage(WavefrontAlignmentThreadLocalStorage &&) = default;

    //-------------------------------------------------------------------------
    // Destructor.

    ~WavefrontAlignmentThreadLocalStorage() = default;

    //-------------------------------------------------------------------------
    // Member Functions.

    WavefrontAlignmentThreadLocalStorage& operator=(WavefrontAlignmentThreadLocalStorage const &) = default;
    WavefrontAlignmentThreadLocalStorage& operator=(WavefrontAlignmentThreadLocalStorage &&) = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TConfig>
inline auto &
intermediate(WavefrontAlignmentThreadLocalStorage<TConfig> & me,
             size_t const alignId)
{
    SEQAN_ASSERT_LT(alignId, me._mMultiAlignmentThreadLocal.size());
    return std::get<typename TConfig::TIntermediate>(me._mMultiAlignmentThreadLocal[alignId]);
}

template <typename TConfig>
inline auto &
cache(WavefrontAlignmentThreadLocalStorage<TConfig> & me,
      size_t const alignId)
{
    SEQAN_ASSERT_LT(alignId, me._mMultiAlignmentThreadLocal.size());
    return std::get<typename TConfig::TCache>(me._mMultiAlignmentThreadLocal[alignId]);
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_ALIGN_PARALLEL_DP_THREAD_LOCAL_STORAGE_H_
