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

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_HORSPOOL_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_HORSPOOL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TFinder>
struct HorspoolFunctor_;

template <typename TContainer, typename TNeedle, typename TSpec>
struct HorspoolFunctor_<Finder2<TContainer, Pattern<TNeedle, Horspool>, TSpec> >
{

    typedef typename Iterator<TNeedle, Rooted>::Type TNeedleIt;

    Pattern<TNeedle, Horspool>* _pattern;  // We keep a pointer to the pattern to access its data map.
    TNeedleIt _itBegin;
    TNeedleIt _itEnd;

    HorspoolFunctor_() : _pattern(NULL), _itBegin(), _itEnd()
    {}

    HorspoolFunctor_(Pattern<TNeedle, Horspool> & pattern)
    {
        _init(*this, pattern);
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

template <typename TNeedle>
struct PatternSpecificTraversalSpec<Pattern<TNeedle, Horspool> >
{
    typedef TraverserSpec<ContextPositionRight, ContextBeginRight> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FinderFunctor                                        [Horspool]
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec>
struct FinderFunctor<Finder2<THaystack, Pattern<TNeedle, Horspool>, TSpec> >
{
    typedef Finder2<THaystack, Pattern<TNeedle, Horspool>, TSpec> TFinder;
    typedef HorspoolFunctor_<TFinder> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _init()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TNeedle>
inline void
_init(HorspoolFunctor_<TFinder> & horspoolFunctor,
      Pattern<TNeedle, Horspool> & pattern)
{
    _patternInit(pattern);  // Initialize the pattern.
    horspoolFunctor._pattern = &pattern;
    horspoolFunctor._itBegin = begin(needle(pattern), Rooted());
    horspoolFunctor._itEnd = end(needle(pattern), Rooted());
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_HORSPOOL_H_
