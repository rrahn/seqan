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
// Implements the default virtual string tree to traverse multiple sequences
// in parallel. This is a facade combining the variant store with the
// journal set.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_JOURNALED_STRING_TREE_DEFAULT_H_
#define EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_JOURNALED_STRING_TREE_DEFAULT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename T>
struct JournalData{};

template <typename T>
struct VariantData{};

// ----------------------------------------------------------------------------
// Class JournaledStringTree                                  [StringTreeDefault]
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore>
class JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault>
{
public:
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore> TDeltaMap;
    typedef typename JournalData<JournaledStringTree>::Type TJournalData;
    typedef typename Size<TJournalData>::Type TSize;

    // TODO(rmaerker): Maybe no holder.
    Holder<TDeltaMap> _variantData;
    Holder<TJournalData>  _journalData;
    TSize _journalSize;

    bool _emptyJournal;
    unsigned _blockBegin;
    unsigned _blockEnd;

    JournaledStringTree() : _variantData(), _journalData(), _journalSize(0), _emptyJournal(true), _blockBegin(0), _blockEnd(0)
    {}

    template <typename THost>
    JournaledStringTree(THost const & host, TDeltaMap const & varData) : _emptyJournal(true)
    {
        _journalSize = coverageSize(deltaCoverageStore(varData));
        setValue(_variantData, varData);
        setHost(*this, host);
        _blockBegin = 0;
        _blockEnd = length(keys(varData));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
struct Spec<JournaledStringTree<TDeltaMap, TSpec> >
{
    typedef TSpec Type;
};

template <typename TDeltaMap, typename TSpec>
struct Spec<JournaledStringTree<TDeltaMap, TSpec> const> :
    Spec<JournaledStringTree<TDeltaMap, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
struct Position<JournaledStringTree<TDeltaMap, TSpec> > :
    Position<TDeltaMap>{};

template <typename TDeltaMap, typename TSpec>
struct Position<JournaledStringTree<TDeltaMap, TSpec> const> :
    Position<TDeltaMap const>{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
struct Size<JournaledStringTree<TDeltaMap, TSpec> > :
    Size<TDeltaMap>{};

template <typename TDeltaMap, typename TSpec>
struct Size<JournaledStringTree<TDeltaMap, TSpec> const> :
    Size<TDeltaMap const>{};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
struct Host<JournaledStringTree<TDeltaMap, TSpec> >
{
    typedef JournaledStringTree<TDeltaMap, TSpec> TVStringTree_;
    typedef typename JournalData<TVStringTree_>::Type TJournalData_;
    typedef typename Host<TJournalData_>::Type Type;
};

template <typename TDeltaMap, typename TSpec>
struct Host<JournaledStringTree<TDeltaMap, TSpec> const>
{
    typedef JournaledStringTree<TDeltaMap, TSpec> TVStringTree_;
    typedef typename Host<TVStringTree_>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction VariantData
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
struct VariantData<JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> >
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore> Type;
};

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
struct VariantData<JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> const>
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
struct JournalData<JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> >
{
    typedef typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_SNP>::Type TAlphabet_;
    typedef typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_DEL>::Type TSize_;

    typedef String<TAlphabet_, Journaled<Alloc<>, SortedArray> > TString_;
    typedef StringSet<TString_, Owner<JournaledSet> > Type;
};

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
struct JournalData<JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> const>
{
    typedef typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_SNP>::Type TAlphabet_;
    typedef typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_DEL>::Type TSize_;

    typedef String<TAlphabet_, Journaled<Alloc<>, SortedArray> > TString_;
    typedef StringSet<TString_, Owner<JournaledSet> > const Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename THost>
inline void
setHost(JournaledStringTree<TDeltaMap, TSpec> & stringTree,
        THost const & newHost)
{
    setHost(journalData(stringTree), newHost);
}

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
inline typename Host<JournaledStringTree<TDeltaMap, TSpec> >::Type &
host(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
{
    return host(journalData(stringTree));
}

template <typename TDeltaMap, typename TSpec>
inline typename Host<JournaledStringTree<TDeltaMap, TSpec> const>::Type &
host(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
{
    return host(journalData(stringTree));
}

// ----------------------------------------------------------------------------
// Function setDeltaMap()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore>
void setDeltaMap(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault> & stringTree,
                 DeltaMap<TDeltaStore, TDeltaCoverageStore> const & deltaMap)
{
    setValue(stringTree._variantData, deltaMap);
    stringTree._blockBegin = 0;
    stringTree._blockEnd = length(keys(deltaMap));
}

// ----------------------------------------------------------------------------
// Function createJournal()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TParallelTag>
void createJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault> & stringTree,
                   Tag<TParallelTag> tag = Serial())
{
    SEQAN_ASSERT_NOT(empty(host(stringTree)));  // The host must be set before.

    adaptTo(value(stringTree._journalData), value(stringTree._variantData), stringTree._blockBegin,
            stringTree._blockEnd, tag);
    stringTree._emptyJournal = false;
}

// ----------------------------------------------------------------------------
// Function requireJournal()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TParallelTag>
void requireJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault> & stringTree,
                    Tag<TParallelTag> tag)
{
    SEQAN_ASSERT_NOT(empty(host(stringTree)));  // The host must be set before.
    if (stringTree._emptyJournal)
        createJournal(stringTree, tag);
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
void requireJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree)
{
    requireJournal(stringTree, Serial());
}

// ----------------------------------------------------------------------------
// Function reinitJournal()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TParallelTag>
void reinitJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree,
                   Tag<TParallelTag> tag = Serial())
{
    clear(journalData(stringTree).strings);
    clear(journalData(stringTree).limits);
    journalData(stringTree).limitsValid = true;
    createJournal(stringTree, tag);
}

// ----------------------------------------------------------------------------
// Function setBlockBegin()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TPosition>
void setBlockBegin(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree,
                   TPosition const & newBegin)
{
    stringTree._blockBegin = newBegin;
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TPosition>
void setBlockEnd(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree,
                 TPosition const & newEnd)
{
    stringTree._blockEnd = newEnd;
}

// ----------------------------------------------------------------------------
// Function variantData()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
inline typename VariantData<JournaledStringTree<TDeltaMap, TSpec> >::Type &
variantData(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
{
    return value(stringTree._variantData);
}

template <typename TDeltaMap, typename TSpec>
inline typename VariantData<JournaledStringTree<TDeltaMap, TSpec> const>::Type &
variantData(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
{
    return value(stringTree._variantData);
}

// ----------------------------------------------------------------------------
// Function journalData()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
inline typename JournalData<JournaledStringTree<TDeltaMap, TSpec> >::Type &
journalData(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
{
    return value(stringTree._journalData);
}

template <typename TDeltaMap, typename TSpec>
inline typename JournalData<JournaledStringTree<TDeltaMap, TSpec> const>::Type &
journalData(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
{
    return value(stringTree._journalData);
}

// ----------------------------------------------------------------------------
// Function keys()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
inline typename Keys<typename VariantData<JournaledStringTree<TDeltaMap, TSpec> >::Type>::Type &
keys(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
{
    return keys(variantData(stringTree));
}

template <typename TDeltaMap, typename TSpec>
inline typename Keys<typename VariantData<JournaledStringTree<TDeltaMap, TSpec> const>::Type >::Type &
keys(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
{
    return keys(variantData(stringTree));
}

// ----------------------------------------------------------------------------
// Function proxyId()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TPosition>
inline typename Position<JournaledStringTree<TDeltaMap, TSpec> const>::Type
proxyId(JournaledStringTree<TDeltaMap, TSpec> const & stringTree,
        TPosition  pos)
{
    SEQAN_ASSERT_GEQ(pos, static_cast<TPosition>(0));
    SEQAN_ASSERT_LT(pos, length(keys(variantData(stringTree))));

    return bitScanForward(mappedCoverage(variantData(stringTree), pos));
}

// ----------------------------------------------------------------------------
// Function checkProxyId()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TPosition>
inline bool
checkProxyId(JournaledStringTree<TDeltaMap, TSpec> const & stringTree,
             TPosition const & proxyId)
{
    return (proxyId > 0) && (proxyId < lenght(value(stringTree._journalData)));
}

// ----------------------------------------------------------------------------
// Function getProxy()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TPosition>
inline typename Value<typename JournalData<JournaledStringTree<TDeltaMap, TSpec> >::Type >::Type &
getProxy(JournaledStringTree<TDeltaMap, TSpec> & stringTree,
         TPosition const & proxyId)
{
    SEQAN_ASSERT_GEQ(proxyId, static_cast<TPosition>(0));
    SEQAN_ASSERT_LT(proxyId, static_cast<TPosition>(length(value(stringTree._journalData))));

    return value(value(stringTree._journalData), proxyId);
}

template <typename TDeltaMap, typename TSpec, typename TPosition>
inline typename Value<typename JournalData<JournaledStringTree<TDeltaMap, TSpec> const>::Type >::Type &
getProxy(JournaledStringTree<TDeltaMap, TSpec> const & stringTree,
         TPosition const & proxyId)
{
    SEQAN_ASSERT_GEQ(proxyId, static_cast<TPosition>(0));
    SEQAN_ASSERT_LT(proxyId, static_cast<TPosition>(length(value(stringTree._journalData))));

    return value(value(stringTree._journalData), proxyId);
}

// ----------------------------------------------------------------------------
// Function journalSize()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
inline typename Size<JournaledStringTree<TDeltaMap, TSpec> const>::Type
journalSize(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
{
    return stringTree._journalSize;
}

// ----------------------------------------------------------------------------
// Function reset()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec>
inline void
shrink(JournaledStringTree<TDeltaMap, TSpec> const & /*stringTree*/)
{
    // Nothing to do.
}

// ----------------------------------------------------------------------------
// Function updateProxy()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TPosition, typename TPosition2>
inline void
updateProxy(JournaledStringTree<TDeltaMap, TSpec> & /*stringTree*/,
            TPosition const & /*proxyId*/,
            TPosition2 const & /*variantId*/)
{

}   // Nothing to do.

// ----------------------------------------------------------------------------
// Function copyProxy()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TJournal, typename TPosition, typename TPosition2>
inline void
copyProxy(JournaledStringTree<TDeltaMap, TSpec> & /*stringTree*/,
          TJournal const & /*source*/,
          TPosition const & /*newProxyId*/,
          TPosition2 const & /*variantId*/)
{
    // Nothing to do.
}

// ----------------------------------------------------------------------------
// Function createProxy()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TPosition>
inline void
createProxy(JournaledStringTree<TDeltaMap, TSpec> & /*stringTree*/,
            TPosition const & /*newProxyId*/)
{
    // Nothing to do.
}

template <typename TDeltaMap, typename TSpec, typename TPosition, typename TPosition2>
inline void
createProxy(JournaledStringTree<TDeltaMap, TSpec> & /*stringTree*/,
            TPosition const & /*newProxyId*/,
            TPosition2 const & /*variantId*/)
{
    // Nothing to do.
}

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_JOURNALED_STRING_TREE_DEFAULT_H_
