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

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_FILTER_HITS_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_FILTER_HITS_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSize = __int64>
class PigeonholeHits
{
public:
    typedef typename Value<PigeonholeHits>::Type     TValue;
    typedef String<TValue>                           TData;
    typedef typename Iterator<TData, Standard>::Type TIterator;

    String<TValue>  _data;
    TIterator       currIter;
    TIterator       endIter;

    PigeonholeHits() : currIter(), endIter()
    {}
};

// ----------------------------------------------------------------------------
// Class PigeonholeHit
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec = void>
struct PigeonholeHit2
{
    TSize ndlSeqPos;
    TSize ndlSeqId;
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TSize>
struct Value<PigeonholeHits<TSize> >
{
    typedef PigeonholeHit2<TSize, void> Type;
};

template <typename TSize>
struct Value<PigeonholeHits<TSize> const>
{
    typedef PigeonholeHit2<TSize, void> Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TSize, typename TExpand>
void
appendValue(PigeonholeHits<TSize> & hits,
            typename Value<PigeonholeHits<TSize> >::Type val,
            Tag<TExpand> const & expand)
{
    if (end(hits._data, Standard()) == hits.endIter)
        appendValue(hits._data, val, expand);
    else
        *hits.endIter = val;
    ++hits.endIter;
}

template <typename TSize>
inline bool
hasNext(PigeonholeHits<TSize> const & hits)
{
    if (empty(hits._data))
        return false;
    return hits.currIter != hits.endIter;
}

template <typename TSize>
inline typename Value<PigeonholeHits<TSize> >::Type &
getNext(PigeonholeHits<TSize> & hits)
{
    return *(hits.currIter++);
}

template <typename TSize>
inline typename Value<PigeonholeHits<TSize> const>::Type &
getNext(PigeonholeHits<TSize> const & hits)
{
    return *(hits.currIter++);
}

template <typename TSize>
void
clear(PigeonholeHits<TSize> & hits)
{
    if (empty(hits._data))
        return;
    hits.currIter = begin(hits._data, Standard());
    hits.endIter = hits.currIter;
}

template <typename TSize>
void
reinit(PigeonholeHits<TSize> & hits)
{
    clear(hits._data);
    clear(hits.currIter);
    clear(hits.endIter);
}

template <typename TSize>
bool empty(PigeonholeHits<TSize> const & hits)
{
    return !hasNext(hits);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_FILTER_HITS_H_
