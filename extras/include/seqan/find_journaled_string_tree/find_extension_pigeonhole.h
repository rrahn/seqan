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
// Implements a pigeonhole filter for traversing multiple reference
// sequences.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_PIGEONHOLE_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_PIGEONHOLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TFinder, typename TSpec>
class FinderExtensionPointState;

template <typename TFinder, typename TSpec>
class FinderExtensionPointState<TFinder, Pigeonhole<TSpec> >
{
public:
    typedef typename GetPattern<TFinder>::Type TPattern;
    typedef typename TPattern::TShape TShape;

    TShape shape;

    FinderExtensionPointState()
    {}
};

// ----------------------------------------------------------------------------
// Class Pattern
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
class FinderExtensionPoint<TFinder, Pigeonhole<TSpec> >
{
public:

    typedef FinderState<Pigeonhole<TSpec> >                         TFinderState;
    typedef typename GetPattern<TFinder>::Type                      TPattern;
    typedef typename Host<TPattern>::Type                           TIndex;

    TPattern *      _patternPtr;
    TFinderState *  _statePtr;
    __int64         _qGramPrefix;
    bool            _firstHash;

    FinderExtensionPoint()
    {}

    template <typename TResult, typename THystkIt, typename TLogical>
    inline void
    operator()(TResult & res, THystkIt haystackIt)
    {
        typedef typename Fibre<TIndex, QGramSA>::Type const TSA;
        typedef typename Iterator<TSA, Standard>::Type      TSAIter;
        typedef typename Value<TFinderState>::Type          THit;
        typedef typename Position<TFinder>::Type            TPosition;
        typedef typename TPattern::TShape                   TShape;
        typedef typename Value<TShape>::Type                THash;

        haystackIt += _qGramPrefix;  // Set to actual begin of qGram.
        TIndex const &index = host(*_patternPtr);
        THash hashVal;
        if (_firstHash)
        {
            hashVal = hash(_patternPtr->shape, haystackIt);
            _firstHash = false;
        }
        else
        {
            hashVal = hashNext(_patternPtr->shape, haystackIt);
        }

        // all previous matches reported -> search new ones
        clear(*_statePtr);

        TSAIter saBegin = begin(indexSA(index), Standard());
        Pair<unsigned> ndlPos;
        THit hit;

        unsigned bktNo = getBucket(index.bucketMap, hashVal);
        TSAIter occ = saBegin + indexDir(index)[bktNo];
        TSAIter occEnd = saBegin + indexDir(index)[bktNo + 1];

        for(; occ != occEnd; ++occ)
        {
            posLocalize(ndlPos, *occ, stringSetLimits(index));
            hit.ndlSeqPos = getSeqOffset(ndlPos);     // relative needle position.
            hit.ndlSeqNo = getSeqNo(ndlPos);                        // needle seq. number

            if (Pigeonhole<TSpec>::ONE_PER_DIAGONAL && !IsBranchIteator<THystkIt>::VALUE)
            {
                __int64 diag = static_cast<__int64>(position(haystackIt)) - hit.hstkPos + (__int64)_patternPtr->finderPosOffset;
                if (_patternPtr->lastSeedDiag[hit.ndlSeqNo] == diag)
                    continue;
                _patternPtr->lastSeedDiag[hit.ndlSeqNo] = diag;
            }
            hit.ndlSeqLength = sequenceLength(hit.ndlSeqNo, host(*_patternPtr));
//            unsigned ndlLength = sequenceLength(hit.ndlSeqNo, host(_pattern));

//            if (Pigeonhole<TSpec>::HAMMING_ONLY != 0)
//            {
//                hit.bucketWidth = ndlLength;
//            }
//            else
//            {
//                unsigned indels = (unsigned)floor(_pattern._currentErrorRate * ndlLength);
//                hit.bucketWidth = ndlLength + (indels << 1);
//                hit.hstkPos -= indels;
//            }
            appendValue(*_statePtr, hit);
        }

        // biggest position encountered.
        _patternPtr->finderPosNextOffset = _max(_patternPtr->finderPosNextOffset, position(haystackIt) + _patternPtr->maxSeqLen);

        res.i2 = !empty(*_statePtr);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                           [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIndex, typename TFilterSpec, typename TSpec>
struct RegisteredExtensionPoint<Finder2<TContainer, Pattern<TIndex, Pigeonhole<TFilterSpec> >, Jst<TSpec> > >
{
    typedef Finder2<TContainer, Pattern<TIndex, Pigeonhole<TFilterSpec> >, Jst<TSpec> > TFinder_;
    typedef FinderExtensionPoint<TFinder_, Pigeonhole<TFilterSpec> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
struct GetState<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > >
{
    typedef typename GetPattern<TFinder>::Type TPattern_;
    typedef typename TPattern_::TShape Type;
};

template <typename TFinder, typename TSpec>
struct GetState<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > const>
{
    typedef typename GetPattern<TFinder>::Type TPattern_;
    typedef typename TPattern_::TShape const Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                 [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
struct RequireFullContext<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > > : False{};

//// ----------------------------------------------------------------------------
//// Metafunction ContextIteratorPosition                            [Pigeonhole]
//// ----------------------------------------------------------------------------
//
//template <typename TFinder, typename TSpec>
//struct ContextIteratorPosition<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > >
//{
//    typedef ContextPositionMiddle Type;
//};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                               [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
inline typename GetState<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > >::Type &
getState(FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > & extensionFunctor)
{
    return extensionFunctor._patternPtr->shape;
}

template <typename TFinder, typename TSpec>
inline typename GetState<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > const>::Type &
getState(FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > const & extensionFunctor)
{
    return extensionFunctor._patternPtr->shape;
}

// ----------------------------------------------------------------------------
// Function setState()                                             [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
inline void
setState(FinderExtensionPoint<TFinder, Pigeonhole<TSpec> >  & extensionFunctor,
         typename GetState<FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > >::Type const & state)
{
    extensionFunctor._patternPtr->shape = state;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
inline Pair<unsigned>
init(FinderExtensionPoint<TFinder, Pigeonhole<TSpec> > & pigeonholeFunctor,
     typename GetPattern<TFinder>::Type & pattern,
     FinderState<Pigeonhole<TSpec> > & state,
     double const & errorRate)
{
    typedef typename GetPattern<TFinder>::Type TPattern;
    typedef typename Host<TPattern>::Type      TIndex;

    _patternInit(pattern, errorRate);

    pigeonholeFunctor._pattern = pattern;
    pigeonholeFunctor._finderState = &state;

    reinit(state);

    unsigned errors = (unsigned) floor(errorRate * pattern.maxSeqLength);
    unsigned qGram = length(indexShape(host(pattern)));
    pigeonholeFunctor._qGramPrefix = pattern.maxSeqLength + errors - (qGram - 1);
    pigeonholeFunctor._firstHash = true;
    return Pair<unsigned>(((pattern.maxSeqLength + errors) << 1) - qGram, (pattern.maxSeqLength + errors) - qGram);
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_PIGEONHOLE_H_
