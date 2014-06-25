// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Implements simple online serarch.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SIMPLE_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SIMPLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FinderFinderExtensionPoint
// ----------------------------------------------------------------------------

template <typename TPattern_>
class FinderExtensionPoint<TPattern_, Simple> : public FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_> TSuper;
    typedef typename Host<TPattern_>::Type THost;
    typedef typename Iterator<THost, Standard>::Type TNeedleIt;

    TNeedleIt _itBegin;
    TNeedleIt _itEnd;

    FinderExtensionPoint(TPattern_ & pattern) : TSuper(pattern)
    {
        init(*this);
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt haystackIt)
    {
        TNeedleIt ndlIt = _itBegin;
        for (; ndlIt != _itEnd; ++ndlIt, ++haystackIt)
            if (*ndlIt != getValue(haystackIt))
                return;
        res.i1 = true;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                               [Simple]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, Simple> >
{
    typedef FinderExtensionPoint<Pattern<TNeedle, Simple>, Simple> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TPattern>
inline void
init(FinderExtensionPoint<TPattern, Simple> & simpleFunctor)
{
    if (isInit(simpleFunctor))
        return;

    simpleFunctor._itBegin = begin(needle(simpleFunctor._pattern), Standard());
    simpleFunctor._itEnd = end(needle(simpleFunctor._pattern), Standard());
    setInit(simpleFunctor);
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SIMPLE_H_
