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

template <typename TPattern_>
class FinderExtensionPoint<TPattern_,  Horspool> : public FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_> TSuper;
    typedef typename Needle<TPattern_>::Type TNeedle;
    typedef typename Iterator<TNeedle, Rooted>::Type TNeedleIt;

    TNeedleIt _itBegin;
    TNeedleIt _itEnd;

    FinderExtensionPoint(TPattern & pattern) : TSuper(pattern);
    {
        init(*this);
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt const & haystackIt)
    {
        THystkIt hystkIt = haystackIt;
        TNeedleIt ndlIt = _itEnd;
        res.i2 = this->_pattern.data_map[ordValue(getValue(hystkIt))];
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

template <typename TPattern_>
struct ContextIteratorPosition<FinderExtensionPoint<TPattern_, Horspool> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                             [Horspool]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, Horspool> >
{
    typedef FinderExtensionPoint<Pattern<TNeedle, Horspool>, Horspool> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline void
init(FinderExtensionPoint<TPattern_, Horspool> & extension)
{
    if (!isInit(extension))
        return;

    _patternInit(getPattern(extension));  // Initialize the pattern.
    extension._itBegin = begin(host(getPattern(extension)), Rooted());
    extension._itEnd = end(host(getPattern(extension)), Rooted());
    setInit(extension);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_HORSPOOL_H_
