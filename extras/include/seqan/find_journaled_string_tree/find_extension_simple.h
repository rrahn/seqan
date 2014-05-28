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

template <typename TContainer, typename TNeedle, typename TSpec>
class FinderExtensionPoint<Finder_<TContainer, Pattern<TNeedle, Simple>, TSpec>, Simple>
{
public:

    typedef typename Iterator<TNeedle, Standard>::Type TNeedleIt;

    TNeedleIt _itBegin;
    TNeedleIt _itEnd;

    FinderExtensionPoint()
    {}

    FinderExtensionPoint(Pattern<TNeedle, Simple> & pattern)
    {
        init(*this, pattern);
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

template <typename TContainer, typename TNeedle, typename TSpec>
struct RegisteredExtensionPoint<Finder_<TContainer, Pattern<TNeedle, Simple>, Jst<TSpec> > >
{
    typedef Finder_<TContainer, Pattern<TNeedle, Simple>, Jst<TSpec> > TFinder_;
    typedef FinderExtensionPoint<TFinder_, Simple> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
init(FinderExtensionPoint<TFinder, Simple> & simpleFunctor,
     typename GetPattern<TFinder>::Type & pattern)
{
    simpleFunctor._itBegin = begin(needle(pattern), Standard());
    simpleFunctor._itEnd = end(needle(pattern), Standard());
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SIMPLE_H_
