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
// Implements the verifier for filtered matches.
// ==========================================================================

#ifndef EXTRAS_APPS_JST_MAPPER_JST_MAPPER_VERIFIER_H_
#define EXTRAS_APPS_JST_MAPPER_JST_MAPPER_VERIFIER_H_

#include "matches.h"
#include "extender_left.h"
#include "extender_right.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TDistanceModel>
class ExtenderState
{
public:
    unsigned                 contigId;
    unsigned                 maxErrorsPerRead;
    unsigned                 minErrorsPerRead;
    unsigned                 seedLength;
    unsigned                 errors;
    Match<MultiGenomeMatch>  currMatch;  // The current match.

    ExtenderState(unsigned contigId, unsigned maxErrorsPerRead, unsigned minErrorsPerRead,
                  unsigned seedLength) :
        contigId(contigId),
        maxErrorsPerRead(maxErrorsPerRead),
        minErrorsPerRead(minErrorsPerRead),
        seedLength(seedLength),
        errors(0)
    {}
};

template <typename TFragmentStore, typename TDelegate, typename TDistanceModel>
class Extender
{
public:
    typedef ExtenderState<TDistanceModel> TExtenderState;

    TFragmentStore &               store;
    TDelegate &                    delegate;
    ExtenderState<TDistanceModel>  extenderState;
    ExtenderLeft<TExtenderState>   extenderLeft;
    ExtenderRight<TExtenderState>  extenderRight;
    bool                           disabled;


    Extender(TFragmentStore & fragStore, TDelegate & delegate, unsigned contigId, unsigned maxErrorsPerRead,
             unsigned minErrorsPerRead, unsigned seedLength, bool disabled) :
                 store(fragStore),
                 delegate(delegate),
                 extenderState(contigId, maxErrorsPerRead, minErrorsPerRead, seedLength),
                 extenderLeft(extenderState),
                 extenderRight(extenderState),
                 disabled(disabled)
    {}

    template <typename TTraversalState, typename TFilterState>
    inline void operator()(TTraversalState const & traversalState, TFilterState const & filterState)
    {
        typedef typename Value<TFilterState>::Type TFilterHit;


        while(hasNext(filterState))
        {
            TFilterHit hit = getNext(filterState);

            // Extract the match -> readNo, relReadPos
            if (isMaster(traversalState))
            {
                if (!_extendHit(*this, hit.ndlSeqPos, hit.ndlSeqId,
                                clippedContextBeginPosition(traversalState, StateTraverseMaster()),
                                contextView(traversalState, StateTraverseMaster()),
                                traversalState), 0)
                    continue;
            }
            else
            {
                if (!_extendHit(*this, hit.ndlSeqPos, hit.ndlSeqId,
                                clippedContextBeginPosition(traversalState, StateTraverseBranch()),
                                contextView(traversalState, StateTraverseBranch()),
                                traversalState), 0)
                    continue;
            }
            // unsigned errorCount = verifyPrefix(contextBegin(finder) + relNedlPos);
            while (hasNext(extenderRight.matches))
            {

                onMatch(extenderState.delegate, )
            }
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TExtender, typename TReadPos, typename TReadId, typename TContextPos, typename TContextView,
          typename TTraversalState>
inline bool _extendHit(TExtender & extender,
                       TReadPos readPos,
                       TReadId readId,
                       TContextPos contigBeginPos,
                       TContextView contigView,
                       TTraversalState const & traverserState,
                       unsigned seedErrors)
{
    typedef typename TExtender::TFragmentStore          TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore      TReadStore;
    typedef typename Value<TReadStore>::Type            TReadSeq;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;

    if (extender.disabled)
        return false;

    TReadSeq & read = extender.store.readSeqStore[readId];
    TReadSeqSize readLength = length(read);

    // Extend left.
    if (!extend(extender.extenderLeft, readPos, readId, contigBeginPos, contigView))
        return false;

    // This removes some duplicates. -> How? If there are 0 errors in the prefix step?
    if (extender.extenderState.errors - seedErrors < extender.extenderState.minErrorsPerRead)
        return false;

    // Extend right.
    if (!extend(extender.extenderRight, readPos, readId, contigBeginPos, traverserState))
        return false;

    // What happens here ...
    // call the right extension.
    // Right extension forces own traversal with own HammingFinder.
    // Explores all different paths.
    // What is the general begin position than?

    // That is seed begin position - contigBeginPos => prefixOffset
    // Each match contained in the right extender must store it's relative begin position.

    // relBeginPos - seedLength - prefixOffset => globalBegin Position.


//    TContigSeqSize matchEnd = contigBegin + extender.seedLength;
//
//    if (seedBegin + extender.seedLength < readLength)
//    {
//        TContigSeqSize contigRightEnd = extender.contigSizes[contigId];
//        if (contigRightEnd > contigBegin + readLength - seedBegin)
//            contigRightEnd = contigBegin + readLength - seedBegin;
//
//        TContigInfix contigRight(contig, contigBegin + extender.seedLength, contigRightEnd);
//        TReadInfix readRight(read, seedBegin + extender.seedLength, readLength);
//
//        if (!_extend(extender, contigRight, readRight, errors))
//            return false;
//
//        matchEnd = contigRightEnd;
//    }
//
//    // This removes some duplicates.
//    if (errors < extender.minErrorsPerRead)
//        return false;
//
//    bool reverseComplemented = _fixReverseComplemented(extender, readId);
//    onMatch(extender.matchesDelegate, contigId, matchBegin, matchEnd, readId, errors, reverseComplemented);

    return true;
}

}

#endif // EXTRAS_APPS_JST_MAPPER_JST_MAPPER_VERIFIER_H_
