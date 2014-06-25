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
// Implements match class and additional utilities to deal with matches.
// ==========================================================================

#ifndef EXTRAS_APPS_JST_MAPPER_MATCHES_H_
#define EXTRAS_APPS_JST_MAPPER_MATCHES_H_

#include <extras/apps/masai/matches.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct MultiGenomeMatch_;
typedef Tag<MultiGenomeMatch_> MultiGenomeMatch;

template <>
struct Match<MultiGenomeMatch> : public Match<>
{
    typedef Match<> TSuper;

    String<bool, Packed<> > coverage;

    Match() : TSuper(), coverage()
    {}
};

template <typename TMatch>
struct MatchBuffer
{
    String<TMatch> matches;
    unsigned beginPos;
    unsigned endPos;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================



// ----------------------------------------------------------------------------
// Function fill()                                                      [Match]
// ----------------------------------------------------------------------------

template <typename TContigId, typename TContigPos, typename TReadId, typename TErrors, typename TCoverage>
inline void fill(Match<MultiGenomeMatch> & match,
                 TContigId contigId,
                 TContigPos beginPos,
                 TContigPos endPos,
                 TReadId readId,
                 TErrors errors,
                 TCoverage const & coverage,
                 bool reverseComplemented)
{
    match.coverage = coverage;
    fill(match, contigId, beginPos, endPos, readId, errors, reverseComplemented);
}

}

#endif // EXTRAS_APPS_JST_MAPPER_MATCHES_H_
