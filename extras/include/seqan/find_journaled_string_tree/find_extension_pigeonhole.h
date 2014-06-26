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


template <typename TPattern_, typename THits>
struct PigeonholeState
{
public:
    typedef typename TPattern_::TShape TShape;

    THits* _hitsPtr;
    TShape shape;

    PigeonholeState() : _hitsPtr(NULL)
    {}
};

// ----------------------------------------------------------------------------
// Class Pattern
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
class FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > : FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_>             TSuper;
    typedef typename GetState<FinderExtensionPoint>::Type   TState;
    typedef typename Host<TPattern_>::Type                  TIndex;
    typedef typename Size<TIndex>::Type                     TSize;
    typedef String<PigeonholeHit2<__int64> >                THits;

    THits   _hits;
    TState  _state;
    TSize   _nonRepeatBegin;
    TSize   _seedLength;

    template <typename TErrorRate>
    FinderExtensionPoint(TPattern_ & pattern, TErrorRate errorRate) : TSuper(pattern)
    {
        init(*this, errorRate);
    }

    template <typename TResult, typename THystkIt, typename TPosition>
    inline void
    operator()(TResult & res, THystkIt haystackIt, TPosition pos)
    {
        typedef typename Fibre<TIndex, QGramSA>::Type const TSA;
        typedef typename Iterator<TSA, Standard>::Type      TSAIter;
        typedef typename Value<THits>::Type                 THit;
        typedef typename TPattern_::TShape                  TShape;
        typedef typename Value<TShape>::Type                THash;

        haystackIt -= _seedLength;
        TIndex const &index = host(getPattern(*this));
        THash hashVal;
        if (pos == _nonRepeatBegin)  // First hash of new non-repeat region.
            hashVal = hash(_state.shape, haystackIt);
        else  // Next hash.
            hashVal = hashNext(_state.shape, haystackIt);

        // all previous matches reported -> search new ones
        clear(_hits);

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
            hit.ndlSeqId = getSeqNo(ndlPos);          // needle seq. number

            appendValue(_hits, hit);

//            if (Pigeonhole<TSpec>::ONE_PER_DIAGONAL && !IsBranchIteator<THystkIt>::VALUE)
//            {
//                __int64 diag = static_cast<__int64>(position(haystackIt)) - hit.hstkPos + (__int64)_patternPtr->finderPosOffset;
//                if (_patternPtr->lastSeedDiag[hit.ndlSeqNo] == diag)
//                    continue;
//                _patternPtr->lastSeedDiag[hit.ndlSeqNo] = diag;
//            }
//            hit.ndlSeqLength = sequenceLength(hit.ndlSeqNo, host(*_patternPtr));
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
        }

        // biggest position encountered.
        getPattern(*this).finderPosNextOffset = _max(getPattern(*this).finderPosNextOffset, position(haystackIt) + getPattern(*this).maxSeqLen);

        res.i1 = !empty(_hits);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                           [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TIndex, typename TFilterSpec>
struct RegisteredExtensionPoint<Pattern<TIndex, Pigeonhole<TFilterSpec> > >
{
    typedef FinderExtensionPoint<Pattern<TIndex, Pigeonhole<TFilterSpec> >, Pigeonhole<TFilterSpec> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
struct GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > >
{
    typedef typename FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> >::THits THits_;
    typedef PigeonholeState<TPattern_, THits_> Type;
};

template <typename TPattern_, typename TSpec>
struct GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > const>
{
    typedef typename FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> >::THits THits_;
    typedef PigeonholeState<TPattern_, THits_> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                 [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
struct RequireFullContext<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > > : False{};

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition                            [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
struct ContextIteratorPosition<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > >
{
    typedef ContextPositionRight Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                               [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
inline typename GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > >::Type &
getState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > & extensionFunctor)
{
    return extensionFunctor._state;
}

template <typename TPattern_, typename TSpec>
inline typename GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > const>::Type &
getState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > const & extensionFunctor)
{
    return extensionFunctor._state;
}

// ----------------------------------------------------------------------------
// Function setState()                                             [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
inline void
setState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> >  & extensionFunctor,
         typename GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > >::Type const & state)
{
    extensionFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function initState()
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
inline void
initState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > & extensionFunctor)
{
    extensionFunctor._state.shape = getPattern(extensionFunctor).shape;
    extensionFunctor._state._hitsPtr = &extensionFunctor._hits;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
inline void
init(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > & extension,
     double errorRate)
{
    if (isInit(extension))
        return;

    _patternInit(getPattern(extension), errorRate);

    initState(extension);

    unsigned errors = (unsigned) floor(errorRate * getPattern(extension).maxSeqLength);
//    unsigned qGram = length(indexShape(host(getPattern(extension))));
    extension._seedLength = length(getPattern(extension).shape);
    extension._nonRepeatBegin = extension._seedLength;
    setContextSize(extension, getPattern(extension).maxSeqLength + errors);
    setInit(extension);
}

// ----------------------------------------------------------------------------
// Function deliverContext()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIndex, typename TFilterSpec, typename TSpec, typename TDelegate,
          typename TTraverser, typename TTag>
inline typename Size<TTraverser>::Type
deliverContext(Finder_<TContainer, Pattern<TIndex, Pigeonhole<TFilterSpec> >, JstFind<TSpec> > & finder,
               TDelegate & delegateFunctor,
               TTraverser & traverser,
               TTag const & /*traverserState*/)
{
    typedef typename Size<TContainer>::Type TSize;

    Pair<bool, TSize> res(false, 1);
    if (_contextEndPosition(traverser, TTag()) <
        finder._extensionFunctor._nonRepeatBegin + finder._extensionFunctor._seedLength)
        return res;

    finder._extensionFunctor(res, contextIterator(traverser, TTag()), _contextEndPosition(traverser, TTag()));

#ifdef DEBUG_DATA_PARALLEL
    if (res.i1)  // Interrupt: Call the DelegateFunctor.
    {
        _printContext(traverser);
        delegateFunctor(traverser, getState(finder._extensionFunctor));
    }
#else
    if (res.i1)  // Interrupt: Call the DelegateFunctor.
        delegateFunctor(traverser, getState(finder._extensionFunctor));
#endif
    // Return to the traverser and continue.
    return res.i2;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_PIGEONHOLE_H_
