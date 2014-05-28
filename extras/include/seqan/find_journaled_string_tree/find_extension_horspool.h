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
// Implements the pattern state for the horspool algorithm.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_HORSPOOL_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_HORSPOOL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TFinder>
class FinderExtensionPoint<TFinder,  Horspool>
{
public:

    typedef typename GetPattern<TFinder>::Type TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    typedef typename Iterator<TNeedle, Rooted>::Type TNeedleIt;

    Pattern<TNeedle, Horspool>* _pattern;  // We keep a pointer to the pattern to access its data map.
    TNeedleIt _itBegin;
    TNeedleIt _itEnd;

    FinderExtensionPoint() : _pattern(NULL), _itBegin(), _itEnd()
    {}

    FinderExtensionPoint(TPattern & pattern)
    {
        init(*this, pattern);
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt const & haystackIt)
    {
        THystkIt hystkIt = haystackIt;
        TNeedleIt ndlIt = _itEnd;
        res.i2 = _pattern->data_map[ordValue(getValue(hystkIt))];
        while(position(ndlIt) > 0)
        {
            if (*(--ndlIt) != getValue(hystkIt))
                return;
            --hystkIt;
        }
        res.i1 = true;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction PatternSpecificTraversalSpec                         [Horspool]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct ContextIteratorPosition<FinderExtensionPoint<TFinder, Horspool> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                             [Horspool]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec>
struct RegisteredExtensionPoint<Finder_<TContainer, Pattern<TNeedle, Horspool>, Jst<TSpec> > >
{
    typedef Finder_<TContainer, Pattern<TNeedle, Horspool>, Jst<TSpec> > TFinder_;
    typedef FinderExtensionPoint<TFinder_, Horspool> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
init(FinderExtensionPoint<TFinder, Horspool> & horspoolFunctor,
     typename GetPattern<TFinder>::Type & pattern)
{
    _patternInit(pattern);  // Initialize the pattern.
    horspoolFunctor._pattern = &pattern;
    horspoolFunctor._itBegin = begin(needle(pattern), Rooted());
    horspoolFunctor._itEnd = end(needle(pattern), Rooted());
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_HORSPOOL_H_
