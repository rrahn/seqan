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

// ----------------------------------------------------------------------------
// Struct FinderInitializationState                                [Pigeonhole]
// ----------------------------------------------------------------------------

template <typename TSpec>
struct FinderInitializationState<Pigeonhole<TSpec> >
{
    typedef Pair<unsigned, unsigned> THit;

    double errorRate;
    String<THit> _currHits;

};

// ----------------------------------------------------------------------------
// Class Pattern
// ----------------------------------------------------------------------------

// Adaption of the original pigeonhole pattern.
template <typename TIndex, typename TSpec>
class Pattern<TIndex, Jst<Pigeonhole<TSpec> > > : public Pattern<TIndex, Pigeonhole<TSpec> >
{
public:

    typedef typename Position<TIndex>::Type TPosition;
    typedef Pair<TPosition, TPosition> THit;
    typedef String<THit> THitString;

    THitString _hitBatch;  // All current hits.

    Pattern() : Pattern<TIndex, Pigeonhole<TSpec> >()
    {
        clear(_hitBatch);
    }

    Pattern(TIndex &_index): Pattern<TIndex, Pigeonhole<TSpec> >(_index)
    {
        clear(_hitBatch);
    }

    Pattern(TIndex const &_index): Pattern<TIndex, Pigeonhole<TSpec> >(_index)
    {
        clear(_hitBatch);
    }

};

// ----------------------------------------------------------------------------
// Struct DataParallelPigeonholeHit
// ----------------------------------------------------------------------------

template <typename TPosition>
struct DataParallelPigeonholeHit
{
    TPosition _needleId;
    TPosition _ndlOffset;
};


template <typename TFinder>
struct FunctorPigeonholeFilter_;

template <typename THaystack, typename TIndex, typename TSpec>
struct FunctorPigeonholeFilter_<Finder_<THaystack, Pattern<TIndex, Jst<Pigeonhole<TSpec> > >, Jst<Pigeonhole<TSpec> > > >
{
    typedef Pattern<TIndex, Jst<Pigeonhole<TSpec> > > TPattern;
    typedef Finder_<THaystack, TPattern, Jst<Pigeonhole<TSpec> > > TFinder;

    TPattern & _pattern;

    bool _isFirstQGram;

    FunctorPigeonholeFilter_(TPattern & pattern) :
                     _pattern(pattern),
                     _isFirstQGram(true)
    {
        _patternInit(pattern, 0);
    }

    template <typename TSpec>
    FunctorPigeonholeFilter_(TPattern & pattern,
                             FinderInitializationState<Pigeonhole<TSpec> > const & state) :
                     _pattern(pattern),
                     _isFirstQGram(true)
    {
        _patternInit(pattern, state.errorRate);
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt const & haystackIt)
    {
        typedef typename Fibre<TIndex, QGramSA>::Type const TSA;
        typedef typename Iterator<TSA, Standard>::Type      TSAIter;
        typedef typename TFinder::TPigeonholeHit            THit;
        typedef typename Position<TFinder>::Type TPosition;

        TIndex const &index = host(_pattern);

        // all previous matches reported -> search new ones
        clear(_pattern._hitBatch);

        TSAIter saBegin = begin(indexSA(index), Standard());
        Pair<unsigned> ndlPos;
        THit hit;

        register unsigned bktNo = getBucket(index.bucketMap, hash);
        register TSAIter occ = saBegin + indexDir(index)[bktNo];
        register TSAIter occEnd = saBegin + indexDir(index)[bktNo + 1];

        // Check all occurrences for a hit.
        for(; occ != occEnd; ++occ)
        {
            posLocalize(ndlPos, *occ, stringSetLimits(index));
            appendValue(_pattern._hitBatch, Pair<unsigned, unsigned>(getSeqNo(ndlPos), getSeqOffset(ndlPos)));

            // TODO(rmaerker): Delegate to the verifier to check if a region needs to be scanned again.
//            if (Pigeonhole<TSpec>::ONE_PER_DIAGONAL)
//            {
//                // TODO(rmaerker): How do we get to this information.
//                __int64 diag = hit.hstkPos + (__int64)_pattern.finderPosOffset;
//                if (_pattern.lastSeedDiag[hit.ndlSeqNo] == diag)
//                    continue;
//                _pattern.lastSeedDiag[hit.ndlSeqNo] = diag;
//            }
//            unsigned ndlLength = sequenceLength(hit.ndlSeqNo, host(_pattern));
//
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

//        finder.curHit = begin(finder.hits, Standard());
//        finder.endHit = end(finder.hits, Standard());
//
//        return !empty(finder.hits);
        res.i1 = !empty(_pattern._hitBatch);

    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename THaystack, typename TPattern, typename TSpec>
struct FinderFunctor<Finder_<THaystack, TPattern, Jst<Pigeonhole<TSpec> > > >
{
    typedef FunctorPigeonholeFilter_<Finder_<THaystack, TPattern, Jst<Pigeonhole<TSpec> > > > Type;
};

// ----------------------------------------------------------------------------
// Metafunction InitStateForFinder
// ----------------------------------------------------------------------------

template <typename THaystack, typename TPattern, typename TSpec>
struct InitStateForFinder<Finder_<THaystack, TPattern, Jst<Pigeonhole<TSpec> > > >
{
    typedef FinderInitializationState<Pigeonhole<TSpec> > Type;
};

template <typename THaystack, typename TPattern, typename TSpec>
struct InitStateForFinder<Finder_<THaystack, TPattern, Jst<Pigeonhole<TSpec> > > const>  :
    InitStateForFinder<Finder_<THaystack, TPattern, Jst<Pigeonhole<TSpec> > > >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function needInitHash()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TIndex, typename TPatternSpec, typename TSpec, typename TIterator>
inline bool
needInitHash(Finder_<THaystack, Pattern<TIndex, Pigeonhole<TPatternSpec> >, Jst<Pigeonhole<TSpec> > > const & finder)
{
    return finder._patternComp._isFirstQGram;
}

// ----------------------------------------------------------------------------
// Function compute()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TIndex, typename TPatternSpec, typename TSpec, typename TIterator>
inline typename ComputeState< typename GetTraverserForFinder_<
                              Finder_<THaystack, Pattern<TIndex, Pigeonhole<TPatternSpec> >,
                              Jst<Pigeonhole<TSpec> > > >::Type>::Type
compute(Finder_<THaystack, Pattern<TIndex, Pigeonhole<TPatternSpec> >, Jst<Pigeonhole<TSpec> > > & finder,
        TIterator const & iter)
{
    typedef Finder_<THaystack, Pattern<TIndex, Pigeonhole<TPatternSpec> >, Jst<Pigeonhole<TSpec> > > TFinder;
    typedef typename GetTraverserForFinder_<TFinder>::Type  TTraverser;
    typedef typename ComputeState<TTraverser>::Type         TComputeState;
    typedef typename Fibre<TIndex, QGramShape>::Type        TShape;

    TShape &shape = getPattern(finder).shape;
    if (needInitHash(finder))
        hashInit(shape, iter);
    else
        hashNext(shape, iter);

    TComputeState state(false, 1);
    finder._finderFuntor(state, iter);
    return state;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_PIGEONHOLE_H_
