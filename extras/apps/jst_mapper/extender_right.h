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
// Implements the algorithms to extend a seed to the right.
// ==========================================================================

#ifndef EXTRAS_APPS_JST_MAPPER_EXTENDER_RIGHT_H_
#define EXTRAS_APPS_JST_MAPPER_EXTENDER_RIGHT_H_

#include "extender.h"

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TContainer, typename TPattern, typename TState>
class ExtenderRight
{
public:
//    typedef typename Infix<TReadSeq>::Type                          TReadInfix;
//    typedef Pattern<TReadInfix, HammingPrefix>                      TPattern;
    typedef typename GetJstTraverser<TContainer, TPattern>::Type    TJstTraverser;
//    typedef MatchCollector<TPattern>                                TCollector;

    TState&         extenderState;
    TPattern        pattern;
    TJstTraverser   traverser;

    ExtenderRight(TState & extenderState) : extenderState(extenderState)
    {}

};
// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TContainer, typename TPattern, typename TReadSeqPos, typename TRead, typename TDelegate,
          typename TTraverser>
inline void
extend(ExtenderRight<TContainer, TPattern, ExtenderState<HammingDistance> > & extender,
       TReadSeqPos readBeginPos,
       TRead & read,
       TDelegate & delegate,
       TTraverser const & traverser)
{
    typedef typename Host<TPattern>::Type TNeedle;

    // From here we need to invoke another traversal.
    // Transfer options? FullState(), ContinueCurrent()
    transferState(extender.traverser, traverser, TraverseRightExtend());

    // Set pattern host.
    TNeedle readInfix = infix(read, readBeginPos + extender.extenderState.seedLength, length(read));
    setHost(extender.pattern, readInfix);

    // We need to set the pattern -> and with the pattern the allowed number of errors.
    find(extender.traverser, extender.pattern, delegate,
         -(extender.extenderState.maxErrorsPerRead - extender.extenderState.errors), JstFind<FindPrefix>())

}

#endif // EXTRAS_APPS_JST_MAPPER_EXTENDER_RIGHT_H_
