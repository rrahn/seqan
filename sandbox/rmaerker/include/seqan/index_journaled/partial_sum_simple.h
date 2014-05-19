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
// Implements a simple dynamic partial sum data structure, which computes
// the partial sum in linear time.
// ==========================================================================

#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_PARTIAL_SUM_SIMPLE_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_PARTIAL_SUM_SIMPLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TValue, typename TSpec = Simple>
class PartialSumManager;

template <typename TValue>
class PartialSumManager<TValue, Simple>
{
public:
    typedef Pair<TValue, TValue> TTableValue_;
    String<TTableValue_> _partialSumTable;  // Store the offsets per position.

    PartialSumManager() : _partialSumTable()
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TValue>
struct Cargo<PartialSumManager<TValue, Simple> >
{
    typedef Pair<TValue, TValue> Type;
};

template <typename TValue>
struct Cargo<PartialSumManager<TValue, Simple> const>
{
    typedef Pair<TValue, TValue> const Type;
};

template <typename TValue>
struct Iterator<PartialSumManager<TValue, Simple>, Standard>
{
    typedef PartialSumManager<TValue, Simple> TPartialSumMgr_;
    typedef typename Cargo<TPartialSumMgr_>::Type TCargo_;
    typedef typename Iterator<String<TCargo_>, Standard>::Type Type;
};

template <typename TValue>
struct Iterator<PartialSumManager<TValue, Simple> const, Standard>
{
    typedef PartialSumManager<TValue, Simple> TPartialSumMgr_;
    typedef typename Cargo<TPartialSumMgr_>::Type TCargo_;
    typedef typename Iterator<String<TCargo_> const, Standard>::Type Type;
};

template <typename TCargo>
struct SimplePartialSumLessThanComparator_
{
    SimplePartialSumLessThanComparator_()
    {}

    bool operator()(TCargo const & targetComp, TCargo const & sourceComp)
    {
        return targetComp.i1 < sourceComp.i1;
    }
};

template <typename TCargo>
struct SimplePartialSumAccumulator_
{

    SimplePartialSumAccumulator_()
    {}

    typename Value<TCargo, 2>::Type
    operator()(typename Value<TCargo, 2>::Type const & val, TCargo const & obj)
    {
        return val + obj.i2;
    }
};


// ============================================================================
// Functions
// ============================================================================

template <typename TValue>
inline typename Iterator<PartialSumManager<TValue, Simple>, Standard>::Type
begin(PartialSumManager<TValue, Simple> & partialSumMgr, Standard const & /*tag*/)
{
    return begin(partialSumMgr._partialSumTable);
}

template <typename TValue>
inline typename Iterator<PartialSumManager<TValue, Simple> const, Standard>::Type
begin(PartialSumManager<TValue, Simple> const & partialSumMgr, Standard const & /*tag*/)
{
    return begin(partialSumMgr._partialSumTable);
}

template <typename TValue>
inline typename Iterator<PartialSumManager<TValue, Simple>, Standard>::Type
end(PartialSumManager<TValue, Simple> & partialSumMgr, Standard const & /*tag*/)
{
    return end(partialSumMgr._partialSumTable);
}

template <typename TValue>
inline typename Iterator<PartialSumManager<TValue, Simple> const, Standard>::Type
end(PartialSumManager<TValue, Simple> const & partialSumMgr, Standard const & /*tag*/)
{
    return end(partialSumMgr._partialSumTable);
}

/*
 * If the element exist, then we update its underlying offset?
 */
template <typename TValue, typename TSpec, typename TPos, typename TOffset>
inline void
insert(PartialSumManager<TValue, TSpec> & partialSumMgr,
       TPos const & pos,
       TOffset const & offset)
{
    typedef PartialSumManager<TValue, TSpec> TPartialSumMgr;
    typedef typename Cargo<TPartialSumMgr>::Type TCargo;
    typedef typename Iterator<TPartialSumMgr, Standard>::Type TItereator;

    if (empty(partialSumMgr._partialSumTable))
    {
        appendValue(partialSumMgr._partialSumTable, TCargo(pos, offset));
        return;
    }

    // Finds upper_bound or element itself.
    TItereator it = std::lower_bound(begin(partialSumMgr, Standard()),
                                     end(partialSumMgr, Standard()), TCargo(pos, 0), SimplePartialSumLessThanComparator_<TCargo>());
    // Element exists. Update its offset.
    if (it->i1 == pos)
        it->i2 += offset;
    else  // Element does not exist. Insert before upper bound
        insertValue(partialSumMgr._partialSumTable, position(it, partialSumMgr._partialSumTable), TCargo(pos, offset));
}

template <typename TValue, typename TPos>
inline TValue
partialSum_(PartialSumManager<TValue, Simple> const & partialSumMgr,
           TPos pos)
{
    typedef PartialSumManager<TValue, Simple> TPartialSumMgr;
    typedef typename Cargo<TPartialSumMgr>::Type TCargo;
    typedef typename Iterator<TPartialSumMgr const, Standard>::Type TItereator;

    // Search for the offset.
    TItereator it = std::upper_bound(begin(partialSumMgr, Standard()),
                                     end(partialSumMgr, Standard()), TCargo(pos, 0), SimplePartialSumLessThanComparator_<TCargo>());
    return std::accumulate(begin(partialSumMgr, Standard()), it, 0, SimplePartialSumAccumulator_<TCargo>());
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(PartialSumManager<TValue, Simple> & partialSumMgr)
{
    clear(partialSumMgr._partialSumTable);
}

}  // namespace seqan

#endif  // SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_PARTIAL_SUM_SIMPLE_H_
