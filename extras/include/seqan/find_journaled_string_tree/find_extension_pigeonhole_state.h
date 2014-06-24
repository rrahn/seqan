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
// Implements filter state used for filter algorithms such as pigeonhole
// or swift filter.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_FILTER_STATE_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_FILTER_STATE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TFilterSpec = void>
class FinderState<Pigeonhole<TFilterSpec> >
{
public:
    typedef typename Value<FinderState>::Type TValue;

    String<TValue>  _data;
    unsigned        currPos;
    unsigned        endPos;

    FinderState() : currPos(0), endPos(0)
    {}
};

// ----------------------------------------------------------------------------
// Class PigeonholeHit
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
class PigeonholeHit2
{
    TSize ndlSeqPos;
    TSize ndlSeqNo;
    TSize ndlSeqLength;
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TSpec>
struct Value<FinderState<Pigeonhole<TSpec> > >
{
    typedef PigeonholeHit2<__int64, TSpec> Type;
};

template <typename TSpec>
struct Value<FinderState<Pigeonhole<TSpec> > const>
{
    typedef PigeonholeHit2<__int64, TSpec> Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TSpec,  typename TExpand>
void appendValue(FinderState<Pigeonhole<TSpec> > & matchState,
                 typename Value<FinderState<Pigeonhole<TSpec> > >::Type val,
                 Tag<TExpand> const & expand)
{
    if (length(matchState._data) == matchState.endPos)
        appendValue(matchState._data, val, expand);
    else
        value(matchState._data, matchState.endPos) = val;
    ++matchState.endPos;
}

template <typename TSpec>
bool hasNext(FinderState<Pigeonhole<TSpec> > const & matchState)
{
    return matchState.currPos < matchState.endPos;
}

template <typename TSpec>
typename Value<FinderState<Pigeonhole<TSpec> > >::Type &
getNext(FinderState<Pigeonhole<TSpec> > & matchState)
{
    return matchState._data[matchState.currPos++];
}

template <typename TSize, typename TSpec>
typename Value<FinderState<Pigeonhole<TSpec> > const >::Type &
getNext(FinderState<Pigeonhole<TSpec> > const & matchState)
{
    return matchState._data[matchState.currPos++];
}

template <typename TSpec>
void clear(FinderState<Pigeonhole<TSpec> > & matchState)
{
    matchState.currPos = 0;
    matchState.endPos = 0;
}

template <typename TSpec>
void reinit(FinderState<Pigeonhole<TSpec> > & matchState)
{
    clear(matchState._data);
    matchState.currPos = 0;
    matchState.endPos = 0;
}

template <typename TSpec>
bool empty(FinderState<Pigeonhole<TSpec> > const & matchState)
{
    return !hasNext(matchState);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_FILTER_STATE_H_
