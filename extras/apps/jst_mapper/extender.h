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

template <typename TDistance>
struct GetVerificationAlgorithm{};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec>
class MatchCollector;

template <typename TNeedle>
class MatchCollector<Pattern<TNeedle, HammingPrefix> >
{
public:
    typedef Pattern<TNeedle, HammingPrefix> TPattern;
    typedef typename PatternState<TPattern>::Type TState;
    typedef Match<MultiGenomeMatch> TMatch;
    typedef String<TMatch> TBuffer;

    TMatch  match;
    TBuffer buffer;

    template <typename TTraverser>
    inline void operator()(TTraverser const & traverser, TState const & state)
    {
        match.beginPos = contextBeginPosition(traverser);
        match.endPosDelta = static_cast<int>(contextEndPosition(traverser)) - match.beginPos;
        match.errors = state.errorCount;
        match.coverage = coverage(traverser);
        appendValue(buffer, match);
    }
};

template <typename TDistanceModel>
class ExtenderState
{
public:
    unsigned                 contigId;
    unsigned                 seedLength;
    int                      maxErrorsPerRead;
    int                      minErrorsPerRead;
    int                      errors;
    Match<MultiGenomeMatch>  currMatch;  // The current match.

    ExtenderState(unsigned contigId, unsigned seedLength, unsigned maxErrorsPerRead, unsigned minErrorsPerRead) :
        contigId(contigId),
        seedLength(seedLength),
        maxErrorsPerRead(maxErrorsPerRead),
        minErrorsPerRead(minErrorsPerRead),
        errors(0)
    {}
};

template <typename TFragmentStore, typename TJst, typename TDelegate, typename TDistanceModel>
class Extender
{
public:
    typedef ExtenderState<TDistanceModel> TExtenderState;
    typedef Infix<TReadSeq>::Type         TReadInfix;
    typedef typename GetVerificationAlgorithm<TDistanceModel>::Type TAlgoSpec;
    typedef Pattern<TReadInfix, TAlgoSpec> TVerifier;
    typedef MatchCollector<TAlgoSpec>      TCollector;

    TFragmentStore &               store;
    TJst &                         jst;
    TDelegate &                    delegate;
    ExtenderState<TDistanceModel>  extenderState;
    ExtenderLeft<TExtenderState>   extenderLeft;

    ExtenderRight<TJst, TVerifier, TExtenderState>  extenderRight;
    TCollector                     collector;
    bool                           disabled;


    Extender(TFragmentStore & fragStore, TJst & jst, TDelegate & delegate, unsigned contigId, unsigned seedLength,
             unsigned maxErrorsPerRead, unsigned minErrorsPerRead, bool disabled) :
                 store(fragStore),
                 jst(jst),
                 delegate(delegate),
                 extenderState(contigId, seedLength, maxErrorsPerRead, minErrorsPerRead),
                 extenderLeft(extenderState),
                 extenderRight(extenderState),
                 disabled(disabled)
    {}

    // I need to know that this is the qGram filter and how i can access the data.
    template <typename TTraversalState, typename TFilterState>
    inline void operator()(TTraversalState const & traversalState, TFilterState const & filterState)
    {
        typedef typename Value<TFilterState>::Type TFilterHit;
        typedef typename TCollector::TBuffer TMatchBuffer;
        typedef typename Iterator<TMatchBuffer, Standard>::Type TBufferIterator;


        while(hasNext(filterState))
        {
            TFilterHit hit = getNext(filterState);
            clear(collector.buffer);  // Clear current buffer to store.
            // Extract the match -> readNo, relReadPos
            if (isMasterState(traversalState))
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

            // Pass matches to delegate.
            TBufferIterator it = begin(collector.buffer, Standard());
            TBufferIterator itEnd = end(collector.buffer, Standard());
            for (; it != itEnd; ++it)
                onMatch(extenderState.delegate, *it);
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <>
struct GetVerificationAlgorithm<HammingDistance>
{
    typedef HammingPrefix Type;
};

template <>
struct GetVerificationAlgorithm<EditDistance>
{
    // TODO(rmaerker): Adapt to Myers.
    typedef HammingPrefix Type;
};

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
                       int seedErrors)
{
    typedef typename TExtender::TFragmentStore          TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore      TReadStore;
    typedef typename Value<TReadStore>::Type            TReadSeq;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;

    typedef typename TExtender::TCollector              TCollector;
    typedef typename TCollector::TBuffer                TBuffer;
    typedef typename Value<TBuffer>::Type               TMatch;
    typedef typename Iterator<TBuffer, Standard>::Type  TBufferIterator;

    if (extender.disabled)
        return false;

    TReadSeq & read = extender.store.readSeqStore[readId];
    TReadSeqSize readLength = length(read);

    // Extend left.
    if (!extend(extender.extenderLeft, readPos, read, contigBeginPos, contigView))
        return false;

    // This removes some duplicates. -> How? If there are 0 errors in the prefix step?
    if (extender.extenderState.errors - seedErrors < extender.extenderState.minErrorsPerRead)
        return false;

    // Extend right.
    // Check if seed is at end of read. -> No extension.
    if (readPos + extender.extenderState.seedLength < length(read))
    {
        extend(extender.extenderRight, readPos, read, extender.collector, traverserState);
        // No suffix matches found.
        if (empty(extender.collector.buffer))
            return false;

        TBufferIterator it = begin(extender.collector.buffer, Standard());
        TBufferIterator itEnd = end(extender.collector.buffer, Standard());
        TContigSeqSize contigPrefixDelta = contigBeginPos - contextEndPosition(traverserState);
        for (; it != itEnd; ++it)
        {
            // Refine the matches.
            *it.beginPos -= contigPrefixDelta;
            *it.coverage &= coverage(traverserState);
            *it.readId = readId;
            *it.contigId = extender.extenderState.contigId;
        }
    }

    // We can just add the one hit.
    // TODO(rmaerker): Add single hit.
    TMatch match;
    fillMatch(match, extender.extenderState.contigId, contigBeginPos, contextEndPosition(traverserState),
              readId, coverage(traverserState), false);/
    return true;
}

}

#endif // EXTRAS_APPS_JST_MAPPER_JST_MAPPER_VERIFIER_H_
