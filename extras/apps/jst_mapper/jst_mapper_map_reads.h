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
// Implements strategy to map the reads.
// ==========================================================================

#ifndef EXTRAS_APPS_JST_MAPPER_JST_MAPPER_MAP_READS_H_
#define EXTRAS_APPS_JST_MAPPER_JST_MAPPER_MAP_READS_H_

#include "jst_mapper.h"
#include "jst_mapper_verifier.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TFragmentStore, typename TContigPos, typename TContigValue, typename TVerifierSpec>
JstMapperResult
mapReads(TFragmentStore & fragStore, DeltaMap<TContigPos, TContigValue> & deltaMap, JstMapperOptions const & options)
{
    typedef typename TFragmentStore::TReadSeqStore                        TReadStore;
    typedef typename TFragmentStore::TContigStore                         TContigStore;
//    typedef typename Value<TContigStore>::Type                      TContigStoreElement;
//    typedef typename TContigStoreElement::TContigSeq                TContigSeq;
//    typedef typename Value<TContigSeq>::Type                        TContigSeqAlphabet;
//    typedef typename Position<TContigSeq>::Type                     TContigSeqPosition;
    typedef Index<TReadStore, IndexQGram<Simple, OpenAddressing> >        TIndex;
    typedef Pattern<TIndex, Pigeonhole<> >                                TPattern;
    typedef DeltaMap<TContigPos, TContigValue>                            TDeltaMap;
    typedef JournaledStringTree<TDeltaMap>                                TJst;

    typedef FinderState<Pigeonhole<> >                                    TFilterState;

    typedef JstMapperVerifier<TFilterState, ResultsWriter, TVerifierSpec> TVerifier;

    ResultsWriter writer(options.output);

    // Step 1) Build the index.
    TIndex index(fragStore.readSeqStore);
    TPattern pattern(index);
    TJst jst(fragStore.contigStore[0].seq, deltaMap);
    TFilterState filterState;
    TVerifier verifier(filterState, writer, 3, options.qGram);

    if (options.jstSize != 0)
        setBlockSize(jst, options.jstSize);

    // Step 2) Trigger the search.
    indexRequire(index);
    find(jst, pattern, filterState, verifier, options.errorRate , Jst<Pigeonhole<> >());

}

#endif // EXTRAS_APPS_JST_MAPPER_JST_MAPPER_MAP_READS_H_
