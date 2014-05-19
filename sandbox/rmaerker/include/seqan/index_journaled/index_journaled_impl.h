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
// Implements the delta index.
// ==========================================================================
#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_IMPL_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_IMPL_H_

#include <math.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TPosition>
struct FindState_
{
    typedef String<TPosition> TPositions;

    TPositions _affectedSAPositions;     // Is a sorted list of positions.

    void operator()(TPosition pos)
    {
        if (empty(_affectedSAPositions))
            appendValue(_affectedSAPositions, pos, Exact());
        else
        {
            insertValue(
                _affectedSAPositions,
                std::upper_bound(begin(_affectedSAPositions, Standard()), end(_affectedSAPositions, Standard()), pos) -
                begin( _affectedSAPositions, Standard()),
                pos);
        }
    }
};

//template <typename TIndex>
//struct IndexWrapper_
//{
//    Holder<TIndex>      ref;
//    typename Fibre<TIndex, EsaSA>::Type invSa;
//
//    IndexWrapper_ ()
//    {}
//
//    IndexWrapper_ (TIndex & host) : ref(host)
//    {}
//
//    IndexWrapper_ (TIndex const & host) : ref(host)
//    {}
//
//    IndexWrapper_(IndexWrapper_ & other) : ref(other.ref), invSa(other.invSa)
//    {}
//
//    IndexWrapper_(IndexWrapper_ const & other) : ref(other.ref), invSa(other.invSa)
//    {}
//
//};

// This is the Delta version of the Esa Index
template < typename TText, typename TSpec, typename TSpecDelta>
class Index<TText, IndexJournaled<IndexEsa<TSpec>, TSpecDelta> >
{
public:
    typedef Index<TText, IndexEsa<TSpec> > THostIndex_;
//    typedef IndexWrapper_<THostIndex_> THostIndex;
    typedef typename SAValue<Index>::Type TSAValue;
    typedef PartialSumManager<TSAValue, TSpecDelta> TPartialSumManager_;
    typedef Holder<typename Fibre<Index, FibreText>::Type>    TTextHolder;   // stores the text of the index.

    // TODO(rmaerker): Wrapper used to get the invSa, needs to be removed later.
    // TODO(rmaerker): We do not really need that.
//    THostIndex                                      ref;    // pointer to the reference index?
    // Here we store journaled tables.

    TTextHolder                                     text;

    // We need a special table to store the offsets here.
    TPartialSumManager_                             psMgr;  // global partial sum data structure.
    typename Fibre<Index, FibreSA>::Type            sa;     // Here we store a journaled string.
    typename Fibre<Index, EsaISA>::Type             invSa;  // Stores the inverse Sa as journaled string.
    typename Fibre<Index, EsaLcp>::Type             lcp;    // longest-common-prefix table


//    typename Fibre<Index, EsaLcpe>::Type            lcpe;       // extended lcp table
//    typename Fibre<Index, EsaChildtab>::Type        childtab;   // child table (tree topology)
//    typename Fibre<Index, EsaBwt>::Type             bwt;        // burrows-wheeler table
//    typename Cargo<Index>::Type                     cargo;      // user-defined cargo

    // Note that tables are created on demand.
    Index() {}

    Index(Index &other):
//        ref(other.ref),
        text(other.text),
        sa(other.sa),
        invSa(other.invSa),
        lcp(other.lcp)
//        lcpe(other.lcpe),
//        childtab(other.childtab),
//        bwt(other.bwt),
//        cargo(other.cargo)
    {}

    Index(Index const &other):
//        ref(other.ref),
        text(other.text),
        sa(other.sa),
        invSa(other.invSa),
        lcp(other.lcp)
//        lcpe(other.lcpe),
//        childtab(other.childtab),
//        bwt(other.bwt),
//        cargo(other.cargo)
    {}

    Index(Index<TText, IndexEsa<TSpec> > &_ref, typename Fibre<Index, EsaText>::Type & journalString) :
        text(),
        psMgr(),
        sa(indexSA(_ref)),
        invSa(indexInvSA(_ref)),
        lcp(indexLcp(_ref))
    {
        setValue(text, journalString);  // Set reference to underlying string.
        setPSManager(sa, psMgr);  // psMgr of SA is dependent.
        // TODO(rmaerker): At the moment we create a copy of the psMgr
        createPSManager(invSa, psMgr);  // psMgr of iSA is owner.
        SEQAN_ASSERT_EQ_MSG(&host(value(text)), &indexText(_ref), "Text of reference index does not correspond to text of delta index");  // Check if object ids are equal.
    }

    Index(Index<TText, IndexEsa<TSpec> > const & _ref, typename Fibre<Index, EsaText>::Type & journalString) :
        text(),
        psMgr(),
        sa(indexSA(_ref)),
        invSa(indexInvSA(_ref)),
        lcp(indexLcp(_ref))
    {
        setValue(text, journalString);  // Set reference to underlying string.
        setPSManager(sa, psMgr);
        // TODO(rmaerker): At the moment we create a copy of the psMgr
        createPSManager(invSa, psMgr);
        SEQAN_ASSERT_EQ_MSG(&host(value(text)), &indexText(_ref), "Text of reference index does not correspond to text of delta index");  // Check if object ids are equal.
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetPartialSumManager_
// ----------------------------------------------------------------------------

template <typename TText, typename THostIndex, typename TDeltaSpec>
struct GetPartialSumManager_<Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > >
{
    typedef Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > TDeltaIndex;
    typedef typename TDeltaIndex::TPartialSumManager_ Type;
};

template <typename TText, typename THostIndex, typename TDeltaSpec>
struct GetPartialSumManager_<Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > const>
{
    typedef Index<TText, IndexJournaled<THostIndex, TDeltaSpec> > TDeltaIndex;
    typedef typename TDeltaIndex::TPartialSumManager_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Host<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > >
{
    typedef Index<TText, IndexEsa<TSpec> > Type;
};

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Host<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const>
{
    typedef Index<TText, IndexEsa<TSpec> > const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                                    [Text]
// ----------------------------------------------------------------------------

// TODO(rmaerker): Use default metafunction.
//template <typename TText, typename TSpec, typename TDeltaSpec>
//struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreText>
//{
//    typedef TText Type;
//};
//
//template <typename TText, typename TSpec, typename TDeltaSpec>
//struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreText>
//{
//    typedef TText Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                                 [RawText]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreRawText>
{
    typedef typename Value<TText>::Type TTextValue_;
    typedef typename Spec<TText>::Type TTextSpec_;
    typedef String<TTextValue_, Journaled<TTextSpec_, SortedArray, Alloc<> > > Type;
};

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreRawText>
{
    typedef  Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex_;
    typedef typename Fibre<TDeltaIndex_, FibreRawText>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                            [Suffix Array]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreSA>
{
    typedef  Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex_;
    typedef typename Fibre<typename Host<TDeltaIndex_>::Type, FibreSA>::Type THostSA_;
    typedef typename Spec<THostSA_>::Type THostSASpec_;
    typedef typename SAValue<TDeltaIndex_>::Type TSAValue_;
    typedef typename GetPartialSumManager_<TDeltaIndex_>::Type TPSMgr_;

    typedef JournalIndexFibre<TSAValue_, Journaled<THostSASpec_, SortedArray>,  TPSMgr_> Type;
};

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreSA>
{
    typedef  Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex_;
    typedef typename Fibre<typename Host<TDeltaIndex_>::Type, FibreSA>::Type THostSA_;
    typedef typename Spec<THostSA_>::Type THostSASpec_;
    typedef typename SAValue<TDeltaIndex_>::Type TSAValue_;
    typedef typename GetPartialSumManager_<TDeltaIndex_>::Type TPSMgr_;

    typedef JournalIndexFibre<TSAValue_, Journaled<THostSASpec_, SortedArray>,  TPSMgr_> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                    [Inverse Suffix Array]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreInvSA> :
    Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreSA>{};

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreInvSA> :
    Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreSA>{};

// ----------------------------------------------------------------------------
// Metafunction Fibre                                                     [LCP]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreLcp>
{
    typedef  Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex_;
    typedef typename Fibre<typename Host<TDeltaIndex_>::Type, FibreLcp>::Type TFibreLcp_;
    typedef typename Spec<TFibreLcp_>::Type TFibreLcpSpec_;
    typedef typename Value<TFibreLcp_>::Type TFibreLcpValue_;
    typedef String<TFibreLcpValue_, Journaled<TFibreLcpSpec_, SortedArray> > Type;
};


template <typename TText, typename TSpec, typename TDeltaSpec>
struct Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreLcp>
{
    typedef  Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex_;
    typedef typename Fibre<typename Host<TDeltaIndex_>::Type, FibreLcp>::Type TFibreLcp_;
    typedef typename Spec<TFibreLcp_>::Type TFibreLcpSpec_;
    typedef typename Value<TFibreLcp_>::Type TFibreLcpValue_;
    typedef String<TFibreLcpValue_, Journaled<TFibreLcpSpec_, SortedArray> > const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()                                                     [SA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreSA>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
         FibreSA const & /*tag*/)
{
    return deltaIndex.sa;
}

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreSA>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const & deltaIndex,
         FibreSA const & /*tag*/)
{
    return deltaIndex.sa;
}

// ----------------------------------------------------------------------------
// Function getFibre()                                                    [ISA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreInvSA>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
         FibreInvSA const & /*tag*/)
{
    return deltaIndex.invSa;
}

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreInvSA>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const & deltaIndex,
         FibreInvSA const & /*tag*/)
{
    return deltaIndex.invSa;
}

// ----------------------------------------------------------------------------
// Function getFibre()                                                    [Lcp]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreLcp>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
         EsaLcp const & /*tag*/)
{
    return deltaIndex.lcp;
}

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreLcp>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const & deltaIndex,
         EsaLcp const & /*tag*/)
{
    return deltaIndex.lcp;
}

// ----------------------------------------------------------------------------
// Function getFibre()                                                   [Text]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreText>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
         FibreText const & /*tag*/)
{
    return value(deltaIndex.text);
}

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreText>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const & deltaIndex,
         FibreText const & /*tag*/)
{
    return value(deltaIndex.text);
}

// ----------------------------------------------------------------------------
// Function getFibre()                                                [RawText]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> >, FibreText>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
         FibreRawText const & /*tag*/)
{
    return value(deltaIndex.text);
}

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Fibre<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const, FibreText>::Type &
getFibre(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > const & deltaIndex,
         FibreRawText const & /*tag*/)
{
    return value(deltaIndex.text);
}



// TODO(rmaerker): We need to add this functions later, when we implement the traversals
//template < typename TText, typename TSpec, typename TDeltaSpec>
//void
//_indexRequireTopDownIteration(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > & index)
//{
//    indexRequire(index, EsaSA());
//    indexRequire(index, EsaLcp());
////    indexRequire(index, EsaChildtab());
//}
//
//template < typename TText, typename TSpec >
//void _indexRequireBottomUpIteration(Index<TText, IndexEsa<TSpec> > &index)
//{
//    indexRequire(index, EsaSA());
//    indexRequire(index, EsaLcp());
//}

// ----------------------------------------------------------------------------
// Function indexCreate()                                        [Wrapper call]
// ----------------------------------------------------------------------------

//template <typename TText, typename TSpec>
//inline void
//indexCreate(IndexWrapper_<Index<TText, IndexEsa<TSpec> > > & indexWrapper, FibreInvSA const & /*tag*/)
//{
//    resize(indexWrapper.invSa, length(getFibre(value(indexWrapper.ref), FibreSA())), Exact());
//    for (unsigned i = 0; i < length(getFibre(value(indexWrapper.ref), FibreSA())); ++i)
//        indexWrapper.invSa[saAt(i, value(indexWrapper.ref))] = i;
//}

template <typename TText, typename TSpec, typename TDeltaSpec>
inline typename Host<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > >::Type &
host(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex)
{
    return value(deltaIndex.ref.ref);
}

// ----------------------------------------------------------------------------
// Function indexCreate()                                        [Suffix Array]
// ----------------------------------------------------------------------------

// TODO(rmaerker): indexCreate creates the host index and sets the tables accordingly.
template <typename TText, typename TSpec, typename TDeltaSpec>
inline bool
indexCreate(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex, FibreSA)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex;
    typedef typename DefaultIndexCreator<TDeltaIndex, FibreSA const>::Type TDefConstruction;
    // Here we call the index create delegate on the host of the esa table.
    return indexCreate(deltaIndex, FibreSA(), TDefConstruction());
    // First, create the underlying index and then reinit the journal string.


    indexCreate(deltaIndex, FibreSA(), TDefConstruction());

}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TSpecAlg>
inline bool
indexCreate(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > &index, FibreSA, TSpecAlg const alg)
{
    resize(host(indexSA(index)), length(host(indexRawText(index))), Exact());
    createSuffixArray(host(indexSA(index)), host(indexText(index)), alg);
    reinit(indexSA(index));
    return true;
}

// ----------------------------------------------------------------------------
// Function indexCreate()                                           [LCP Array]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
inline bool
indexCreate(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
            FibreLcp const & /*lcp*/)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex;
    typedef typename DefaultIndexCreator<TDeltaIndex, FibreLcp const>::Type TDefConstruction;

    return indexCreate(deltaIndex, FibreLcp(), TDefConstruction());
}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TSpecAlg>
inline bool
indexCreate(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > &index, FibreLcp, TSpecAlg const alg)
{
    resize(host(indexLcp(index)), length(host(indexRawText(index))), Exact());
    createLcpTable(host(indexLcp(index)), host(indexText(index)), host(indexSA(index)), alg);
    indexLcp(index)._length = length(host(indexRawText(index)));
    reinit(indexLcp(index));
    return true;
}

// ----------------------------------------------------------------------------
// Function indexCreate()                                [Inverse Suffix Array]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec>
void
indexCreate(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
            FibreInvSA const & /*tag*/)
{
    // TODO(rmaerker): Requires the SA and creates it on demand.
    //indexRequire(deltaIndex, FibreSA());
    if (empty(deltaIndex.sa))
        indexCreate(deltaIndex, FibreSA());
    // TODO(rmaerker): Since the underlying index does not know the inv sa we construct it here.
    resize(host(deltaIndex.invSa), length(host(indexText(deltaIndex))), Exact());
    for(unsigned i = 0; i < length(host(deltaIndex.invSa)); ++i)
        host(deltaIndex.invSa)[host(indexSA(deltaIndex))[i]] = i;
    reinit(deltaIndex.invSa);
}

// ----------------------------------------------------------------------------
// Function _findAffectedSuffices()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDeltaSpec, typename TInsertPos, typename TPos, typename TSize>
inline bool
_findAffectedSuffices(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > const & index,
                      TInsertPos insPos,
                      TPos textPos,
                      FindState_<TSize> & findState)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > TIndex;
    typedef typename SAValue<TIndex>::Type TSAValue;

    TSAValue indexPos = invSaAt(textPos, index);
    TSAValue sensitiveIntervalLeft = saAt(indexPos, index) + 1;
    TSAValue sensitiveIntervalRight = saAt(indexPos, index) + _max(lcpAt(indexPos, index), lcpAt(_max(indexPos - 1, 0u), index));
    if (sensitiveIntervalLeft <= insPos && insPos <= sensitiveIntervalRight)
    {
        findState(indexPos);
        return true;
    }
    return false;
}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TSize>
inline void
_deleteAffectedSuffices(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > & index,
                        FindState_<TSize> const & findState)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > TIndex;
    typedef typename FindState_<TSize>::TPositions const TPosArray;
    typedef typename Iterator<TPosArray, Standard>::Type TIter;
    typedef typename SAValue<TIndex>::Type TSAValue;

    // TODO(rmaerker): When we delete the position, we need to know, if it was a position from
    // the insertion buffer or from the original text.
    // If the position came from the insertion buffer, then its value represents already the correct position within
    // the modified index. If not, the value corresponds to the original index, whose values might have a different offset.

    // We delete from right to left, which is easier to handle.
    TIter it = end(findState._affectedSAPositions, Standard());
    do
    {
        --it;
//        std::cerr << "DEBUG: *it == " << *it << std::endl;
//        std::cerr << "DEBUG: length(getFibre(index, FibreLcp())) == " << length(getFibre(index, FibreLcp())) << std::endl;
        // Now we need to update the LCP table accordingly.
        // If the suffix array position is not on the last position.
        if (*it < length(getFibre(index, FibreSA())))
        {
            TSAValue upperLcp = lcpAt(*it - 1, index);
            TSAValue minLcp = _min(upperLcp, lcpAt(*it, index));
            if (minLcp != upperLcp)
                assignLcpValue(indexLcp(index), *it - 1, minLcp);
        }
        // Delete the affected positions in the journaled sequences.
        erase(getFibre(index, FibreLcp()), *it);
        erase(getFibre(index, FibreSA()), *it);
    } while(it != begin(findState._affectedSAPositions, Standard()));
    // TODO(rmaerker): Synchronize ChildTab
        // How to synchronize ChildTab
        // How to synchronize the invSA?
}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TInsertPos, typename TOffset>
void recordOffset(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > & index,
                  TInsertPos const & pos,
                  TOffset const & offset)
{
    insert(index.psMgr, pos, offset);
}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TTextPos, typename TInsertPos,
          typename TInsertText>
inline void
_insertAffectedAndNewSuffices(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > & index,
                              TTextPos textPos,
                              TInsertPos const & insertPos,
                              TInsertText const & insertText)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > TDeltaIndex;
    typedef typename Position<TDeltaIndex>::Type TPos;
    typedef typename Fibre<TDeltaIndex, FibreText>::Type TTextFibre;
    typedef typename Infix<TTextFibre>::Type TInfix;
    typedef Index<TText, IndexSa<> > TSAIndex;

    // Now we have to insert each suffix into the current text.
    // We need to find the occurrence of the substring within the tree.

    // A) We add the offset to the index.
    recordOffset(index, insertPos, length(insertText));

    // B) We scan from left to right over the new and affected suffices and find them in the SA.
    // Old suffixes. Check if they came from the insertion buffer or not.
    // TODO(rmaerker): We now have to update the inverse SA.
    for (TTextPos currPos = textPos; currPos < insertPos; ++currPos)
    {
        LcpAccessor<LcpLowerUpper> lcpAccessor;
        unsigned pos = lowerBoundSA(indexText(index), indexSA(index), suffix(indexText(index), currPos), lcpAccessor);
#ifdef DEBUG_DELTA_INDEX
        std::cout << "Found suffix: " << suffix(indexText(index), currPos) << " at " << pos << std::endl;
        std::cout << "LcpLower: " << lcpAccessor.lcpLower << " LcpUpper: " << lcpAccessor.lcpUpper << std::endl;
        // Update SA and lcp.
        std::cout << "Current text pos: "<< textPos << std::endl;
        if (pos < length(indexSA(index)))
        {
            std::cout << "Current SA value: "<< getValue(indexSA(index), pos) << std::endl;
            std::cout << "Current LCP value: "<< lcpAt(pos, index) << std::endl;
        }
        std::cout << "iSA value:" << invSaAt(textPos, index) << std::endl;
        std::cout << "# Entries to update: " << invSaAt(textPos, index) - pos <<std::endl;
#endif
//        int iSADelta = invSaAt(textPos, index) - pos;
//        // TODO(rmaerker): If iSaDelta < 0 go from position backward otherwise go forward.
//        if (iSADelta < 0)
//        {
//            // TODO(rmaerker): We need a proxy to solve this here.
//            for (unsigned count = 0; count < std::abs(iSADelta); ++count)
//                if (pos >= count)
//                    invSaAt(saAt(pos - count, index), index) += 1;
//            invSaAt(textPos, index) = pos;
//        }
//        else if (iSADelta > 0)
//        {
//            for (unsigned count = 0; count < std::abs(iSADelta); ++count)
//                if (pos >= count)
//                    invSaAt(saAt(pos + count, index), index) += 1;
//            invSaAt(textPos, index) = pos;
//        }
        ++textPos;
        insertValue(indexSA(index), pos, currPos);
        if (pos > 0 && (value(indexLcp(index), pos - 1) != lcpAccessor.lcpLower))
        {
//            std::cout << "Before update: " << value(indexLcp(index), pos - 1) << std::endl;
            assignLcpValue(indexLcp(index), pos - 1, lcpAccessor.lcpLower);
//            value(indexLcp(index), pos - 1) = lcpAccessor.lcpLower;
//            std::cout << "After update: " << value(indexLcp(index), pos - 1) << std::endl;
        }
        insertValue(indexLcp(index), pos, lcpAccessor.lcpUpper);

//        insert(indexSA(index), pos, textPos);
    }

    // New suffixes to add.
    for (TTextPos currPos = insertPos; currPos < insertPos + length(insertText); ++currPos)
        {
            LcpAccessor<LcpLowerUpper> lcpAccessor;
            unsigned pos = lowerBoundSA(indexText(index), indexSA(index), suffix(indexText(index), currPos), lcpAccessor);
#ifdef DEBUG_DELTA_INDEX
            std::cout << "Found suffix: " << suffix(indexText(index), currPos) << " at " << pos << std::endl;
            std::cout << "LcpLower: " << lcpAccessor.lcpLower << " LcpUpper: " << lcpAccessor.lcpUpper << std::endl;
            // Update SA and lcp and iSA.
            std::cout << "Current text pos: "<< textPos << std::endl;
            std::cout << "iSA value:" << invSaAt(textPos, index) << std::endl;
            ++textPos;
#endif
            if (pos < length(indexSA(index)))
            {
                std::cout << "Current SA value: "<< getValue(indexSA(index), pos) << std::endl;
                std::cout << "Current LCP value: "<< lcpAt(pos, index) << std::endl;
            }

            insertValue(indexSA(index), pos, currPos);
            if (pos > 0 && (value(indexLcp(index), pos - 1) != lcpAccessor.lcpLower))
            {
//                std::cout << "Before update: " << value(indexLcp(index), pos - 1) << std::endl;
                assignLcpValue(indexLcp(index), pos - 1, lcpAccessor.lcpLower);
//                std::cout << "After update: " << value(indexLcp(index), pos - 1) << std::endl;
            }
            insertValue(indexLcp(index), pos, lcpAccessor.lcpUpper);
    //        insert(indexSA(index), pos, textPos);
        }


    // Assume we have different tables available.
    // Search costs: O(|P| * log n) with SA
    // Search costs: O(|P| + log n) with SA + LCP // Probably no speed up, but maybe yes, because P is long.
    // Search costs: O(|P|) with SA + LCP + ChildTab
}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TPos, typename TSize>
void
recordDeletion(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
               TPos  deletePos,
               TSize deletionSiyze)
{
    // TODO(rmaerker): Write me!
}

template <typename TText, typename TSpec, typename TDeltaSpec, typename TPos, typename TInsertText>
void
recordInsertion(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex,
                TPos insertPos,
                TInsertText const & insertText)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > TDeltaIndex;
    typedef typename Size<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > >::Type TSize;
    typedef FindState_<TSize> TFinder;
    typedef typename TFinder::TPositions TPosArray;
    typedef typename Iterator<TPosArray, Standard>::Type TIter;
    typedef typename SAValue<TDeltaIndex>::Type TSAValue;

    TFinder findState;
    TPos textPos = insertPos;
    // Step A: Find all entries of SA whose lexicographical order might be affected by the insertion
    bool isAffected = _findAffectedSuffices(deltaIndex, insertPos, --textPos, findState);
    while(isAffected && textPos > 0)
    {
        isAffected = _findAffectedSuffices(deltaIndex, insertPos, --textPos, findState);
    }
    if (textPos > 0)
        ++textPos;

    // TODO(rmaerker): We need to delete the positions in descending order.
    // TODO(rmaerker): When deleting and reinserting the relative positions of the other might change.
    TIter it = end(findState._affectedSAPositions, Standard());
    do
    {
        --it;
        // Now we need to update the LCP table accordingly.
        // If the suffix array position is not on the last position.
        if (*it < length(getFibre(deltaIndex, FibreSA())))
        {
            TSAValue upperLcp = lcpAt(*it - 1, index);
            TSAValue minLcp = _min(upperLcp, lcpAt(*it, index));
            if (minLcp != upperLcp)
                assignLcpValue(indexLcp(deltaIndex), *it - 1, minLcp);
        }
        // Delete the affected positions in the journaled sequences.
        TSAValue saValue = saAt(*it, index);  // We store the suffix to be deleted.

        //now we find the new insertion site of this suffix.
        //updateInvSA

        //updateSA
        erase(getFibre(deltaIndex, FibreLcp()), *it);
        //updateLcp
        erase(getFibre(deltaIndex, FibreSA()), *it);

        //updateLcp

    } while(it != begin(findState._affectedSAPositions, Standard()));
    // TODO(rmaerker): Synchronize ChildTab

    // When we found all affected suffixes:
    // Step a) remove first index from SA
    // Step b) find new SA in index -> wither we find it before another affected suffix > nothing changed. Or after another affected, but its position will be updated anyway
    // Step c) update the invSA for the update position.
    //         Attention: If the affected position is before the insert position nothing changes
                       // But if the new insert position is after the affected position we have deleted the pos from the SA and we need to check the position SA[i-1]
    // Step d) reinsert the new suffix.
    // Step e) update the lcp table.

    // Step B: Delete all entries from the SA.
    // Deletes the affected suffices and updates the lcp table accordingly.
    if (!empty(findState._affectedSAPositions))
        _deleteAffectedSuffices(deltaIndex, findState);
    // TODO(rmaerker): Note, we only need the tables if they are constructed.
    // Step C: Insert affected and new suffices.
    _insertAffectedAndNewSuffices(deltaIndex, textPos, insertPos, insertText);
}

// Note: We assume a sequential order of the operations, otherwise the delta index synchronization becomes to complex
// and probably to expensive.
template <typename TText, typename TSpec, typename TDeltaSpec>
inline void
synchronize(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > & deltaIndex)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec> > TDeltaIndex;
    typedef typename SAValue<TDeltaIndex>::Type TSAValue;
    typedef typename Size<TDeltaIndex>::Type TSize;
    typedef typename Fibre<TDeltaIndex, FibreText>::Type TJournaledString;
    typedef typename JournalType<TJournaledString>::Type TJournalEntries;
    typedef typename Value<TJournalEntries>::Type TJournalEntry;
    typedef typename Position<TJournalEntry>::Type TPos;

    // We require SA, invSA and LCP to be built before.

    // A) Reinit the journaled tables, because we don't know which operations are new.
    reinit(indexSA(deltaIndex));
    reinit(indexLcp(deltaIndex));

    // B) Construct the invSA
    String<TSAValue> invSa;
    resize(invSa, length(host(indexText(deltaIndex))), Exact());
    for(unsigned i = 0; i < length(invSa); ++i)
        invSa[host(indexSA(deltaIndex))[i]] = i;

    // C) Find all affected suffices for all operations and apply changes.
    TJournalEntries & journalEntries = indexText(deltaIndex)._journalEntries;

    TPos lastPhysPos = 0;
    for (unsigned i = 0; i < length(journalEntries._journalNodes); ++i)
    {
        TJournalEntry & entry = journalEntries._journalNodes[i];
        if (entry.segmentSource == SOURCE_ORIGINAL)
        {
            if (lastPhysPos < entry.physicalPosition)
            {  // Record deletion.
                recordDeletion(deltaIndex, lastPhysPos, entry.physicalPosition - lastPhysPos);
            }
            lastPhysPos = entry.physicalPosition + entry.length;
        }
        else if (entry.segmentSource == SOURCE_PATCH)
        {
            std::cout << "The text: " << indexText(deltaIndex) << std::endl;
            std::cout << "The infix: " << infix(indexText(deltaIndex)._insertionBuffer, entry.physicalPosition,
                                                entry.physicalPosition + entry.length) << std::endl;
            recordInsertion(deltaIndex, lastPhysPos,
                            infix(indexText(deltaIndex)._insertionBuffer, entry.physicalPosition,
                                  entry.physicalPosition + entry.length));
        }
    }
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

// TODO(rmaerker): We need a version for StringSets too.
template <typename TText, typename TSpec, typename TDeltaSpec, typename TPos, typename TInsertText>
inline void
insert(Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > & deltaIndex,
       TPos insertPos,
       TInsertText const & insertText)
{
    typedef Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > TDeltaIndex;
    typedef typename Size<Index<TText, IndexJournaled<IndexEsa<TSpec>, TDeltaSpec > > >::Type TSize;

    FindState_<TSize> findState;
    TPos textPos = insertPos;
    // Step a: Find all entries of SA whose lexicographical order might be affected by the insertion
    bool isAffected = _findAffectedSuffices(deltaIndex, insertPos, --textPos, findState);
    while(isAffected)
    {
        isAffected = _findAffectedSuffices(deltaIndex, insertPos, --textPos, findState);
    }

    // Step b: Delete all entries from the SA.
    // Deletes the affected suffices and updates the lcp table accordingly.
    if (!empty(findState._affectedSAPositions))
        _deleteAffectedSuffices(deltaIndex, findState);

    // Step c: Reinsert all new suffices and all old suffices.
    // We have to resize the tables accordingly.
    // TODO(rmaerker): Note, we only need the tables if they are constructed.
    _insertAffectedAndNewSuffices(deltaIndex, ++textPos, insertPos, insertText);

    // Step : Record erase in the underlying text.
}


// ----------------------------------------------------------------------------
// Function getFibre
// ----------------------------------------------------------------------------

}  // namespace seqan

#endif  // SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_IMPL_H_
