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
// Implements a simple dynamic prefix sum data structure, which computes
// the prefix sum in linear time.
// ==========================================================================

#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_INDEX_JOURNALED_INDEX_JOURNALED_BASE_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_INDEX_JOURNALED_INDEX_JOURNALED_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct FibreInvSA_;
typedef Tag<FibreInvSA_> const   FibreInvSA;

typedef FibreInvSA EsaISA;

template <typename THostIndex, typename TSpec = Simple>
struct IndexJournaled;

struct JournaledDeltaIndexFibreIterSpec_;
typedef Tag<JournaledDeltaIndexFibreIterSpec_> JournaledDeltaIndexFibreIterSpec;

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TIndex>
struct GetPartialSumManager_;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function saAt()                                                [Delta Index]
// ----------------------------------------------------------------------------

// TODO(rmaerker): This is not a default implementation.
// TODO(rmaerker): We need to overload these helper functions for the delta index, because the Journaled String Reference Type is a proxy which does not behave as its type.
//template <typename TPosition, typename TText, typename THostIndex, typename TDeltaSpec>
//inline typename Reference<typename Fibre<Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const, FibreSA>::Type>::Type
//saAt_(TPosition pos,
//     Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const & index)
//{
////    static_cast<Nothing>(typename GetValue<typename Fibre<Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const, FibreSA>::Type>::Type());
////    static_cast<Nothing>(typename Fibre<Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const, FibreSA>::Type());
//    return value(getFibre(index, EsaISA()), pos);
//}

// TODO(rmaerker): At the moment this creates a copy.
// Proxy for the inverse sa since it does not exist at the moment.
template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec>, FibreInvSA>::Type
indexInvSA(Index<TText, TSpec> & /*index*/)
{
    return typename Fibre<Index<TText, TSpec>, FibreInvSA>::Type();
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec> const, FibreInvSA>::Type
indexInvSA(Index<TText, TSpec> const & /*index*/)
{
    return typename Fibre<Index<TText, TSpec>, FibreInvSA>::Type();
}

// ----------------------------------------------------------------------------
// Function invSaAt()                                             [Delta Index]
// ----------------------------------------------------------------------------

template <typename TPosition, typename TText, typename THostIndex, typename TDeltaSpec>
inline typename Reference<typename Fibre<Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const, FibreInvSA>::Type >::Type
invSaAt(TPosition pos,
        Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const & index)
{
    return value(getFibre(index, EsaISA()), pos);
}

// ----------------------------------------------------------------------------
// Function assignLcpValue()
// ----------------------------------------------------------------------------

template <typename TLcpFibre, typename TPos, typename TValue>
inline void
assignLcpValue(TLcpFibre & target, TPos const & pos, TValue const & lcpValue)
{
    // Now the target should be a ProxyIterator Type.
    typedef typename Reference<TLcpFibre>::Type TRefValue;

    TRefValue refVal = value(target, pos);
    switch(iter(refVal)._journalEntriesIterator->segmentSource)
    {
        case SOURCE_PATCH:
        {
            value(iter(refVal)._currentInsertionBufferIt) = lcpValue;
            break;
        }
        case SOURCE_ORIGINAL:
        {
            assignValue(iter(refVal), lcpValue);
            break;
        }
        default: SEQAN_ASSERT_FAIL("No valid segment source!");
    }
}

}  // namespace seqan

#endif  // SANDBOX_RMAERKER_INCLUDE_SEQAN_INDEX_JOURNALED_INDEX_JOURNALED_BASE_H_
