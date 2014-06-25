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
// Implements the algorithms to extend a seed to the left.
// ==========================================================================

#ifndef EXTRAS_APPS_JST_MAPPER_EXTENDER_LEFT_H_
#define EXTRAS_APPS_JST_MAPPER_EXTENDER_LEFT_H_

#include "extender.h"

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TState>
class ExtenderLeft;

template <>
class ExtenderLeft<ExtenderState<HammingDistance> >
{
public:
    typedef ExtenderState<HammingDistance> TState;

    TState&      extenderState;

    ExtenderLeft(TState & extenderState) : extenderState(extenderState)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TReadSeqPos, typename TRead, typename TContigPos, typename TContigView>
inline bool
extend(ExtenderLeft<ExtenderState<HammingDistance> > & extender,
       TReadSeqPos readBeginPos,
       TRead & read,
       TContigPos & contigBeginPos,
       TContigView & contigView)
{
    typedef typename Iterator<TContigView, Standard>::Type  TIterContigView;
    typedef typename Segment<TRead const, InfixSegment>     TReadInfix;
    typedef typename Iterator<TReadInfix, Standard>::Type   TIterReadInfix;

    // Seed at the beginning of the read.
    if (readBeginPos == 0)
    {
        contigBeginPos += (end(contigView, Standard()) - extender.extenderState.seedLength) - begin(contigView, Standard());
        return true;
    }

    SEQAN_ASSERT_LEQ(begin(contigView), end(contigView) - extender.extenderState.seedLength - readBeginPos);
    // Clip the end of the contig to the begin of the seed.
    end(contigView, Standard()) -= extender.extenderState.seedLength;

    // Clip the begin of the contig to the left most possible match position of the read.
    if (end(contigView, Standard()) - readBeginPos > begin(contigView, Standard()))
    {
        contigBeginPos += (end(contigView, Standard()) - readBeginPos) - begin(contigView, Standard());
        begin(contigView, Standard()) = end(contigView, Standard()) - readBeginPos;
    }

    TReadInfix readInfix(read, 0, readBeginPos);

    // Check if read is valid.
    if (length(contigView) != length(readInfix))
        return false;

    TIterContigView contigIt = begin(contigView, Standard());
    TIterReadInfix readBegin = begin(readInfix, Standard());
    TIterReadInfix readEnd = end(readInfix, Standard());

    for (TIterReadInfix readIt = readBegin; readIt != readEnd; ++readIt, ++contigIt)
        if (getValue(readIt) != getValue(contigIt))
            if (++extender.extenderState.errors > extender.extenderState.maxErrorsPerRead)
                return false;

    // Prefix matches with at most maxErrorsPerRead!
    return true;
}

}

#endif // EXTRAS_APPS_JST_MAPPER_EXTENDER_LEFT_H_
