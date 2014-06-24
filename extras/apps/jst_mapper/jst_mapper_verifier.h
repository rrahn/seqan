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

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TFilterState, typename TDelegate, typename TSpec>
class JstMapperVerifier
{
public:
    typedef typename Value<TFilterState>::Type TMatchState;

    TMatchState *  statePtr;
    unsigned       maxErrors;
    unsigned       qGram;
    TDelegate      delegate;

    JstMapperVerifier(TMatchState & state, unsigned errors, unsigned qGram) :
        statePtr(&state),
        maxErrors(errors),
        qGram(qGram),
        delegate()
    {}

    template <typename TFinder>
    inline void operator(TFinder const & finder)
    {
        template typename Positions<TFinder>::Type TPositions;

        while(hasNext(*statePtr))
        {
            TMatchState mState = getNext(*statePtr);
            // Extract the match -> readNo, relReadPos

            // unsigned errorCount = verifyPrefix(contextBegin(finder) + relNedlPos);

            if (errorCount < maxErrors)
//                errorCount += verifySuffix(finder, contextBegin(finder) + relNdlPos + qram);
            // Delegate the hits here.
            if (errorCount < maxErrors)
//                delegate(Hit);
        }
        clear(*statePtr);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

inline bool
verifyPrefix(JstMapperVerifier<HammingDistance> & verifier)
{
    return true;
}

inline bool
verifyPrefix(JstMapperVerifier<EditDistance> & verifier)
{
    return true;
}

inline bool
verifySuffix(JstMapperVerifier<HammingDistance> & verifier)
{
    return true;
}

inline bool
verifySuffix(JstMapperVerifier<EditDistance> & verifier)
{
    return true;
}


#endif // EXTRAS_APPS_JST_MAPPER_JST_MAPPER_VERIFIER_H_
