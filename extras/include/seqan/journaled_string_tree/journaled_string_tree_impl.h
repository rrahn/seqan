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
    typedef typename Value<TJournalData>::Type TJournaledString;
    typedef typename Size<TJournalData>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    // TODO(rmaerker): Maybe no holder.
    Holder<TDeltaMap> _variantData;
    Holder<TJournalData>  _journalData;
    String<TSignedSize> _blockVPOffset;  // Virtual position offset for each sequence visited so far.
    String<TSignedSize> _activeBlockVPOffset;  // Virtual position offset for each sequence visited so far.

    static const TSize REQUIRE_FULL_JOURNAL = MaxValue<unsigned>::VALUE;
    TSize _journalSize;
    TSize _blockSize;
    TSize _activeBlock;
    TSize _numBlocks;
    bool _emptyJournal;

    JournaledStringTree() : _variantData(),
                            _journalData(),
                            _journalSize(0),
                            _blockSize(REQUIRE_FULL_JOURNAL),
                            _activeBlock(0),
                            _numBlocks(1),
                            _emptyJournal(true)
    {}

    template <typename THost>
    JournaledStringTree(THost & reference, TDeltaMap const & varData) : _blockSize(REQUIRE_FULL_JOURNAL),
                                                                        _activeBlock(0),
                                                                        _numBlocks(1),
                                                                        _emptyJournal(true)

    {
        _journalSize = coverageSize(deltaCoverageStore(varData));
        setValue(_variantData, varData);
        setHost(*this, reference);
        TJournaledString tmp;
        setHost(tmp, host(*this));
        resize(journalData(*this), _journalSize, tmp);  // Initialize the journaled set.
        resize(_blockVPOffset, _journalSize, 0, Exact());
        resize(_activeBlockVPOffset, _journalSize, 0, Exact());
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

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TMapSpec, typename TSpec>
struct VariantData<JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore, TMapSpec>, TSpec> >
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore, TMapSpec> Type;
};

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TMapSpec, typename TSpec>
struct VariantData<JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore, TMapSpec>, TSpec> const>
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore, TMapSpec> const Type;
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

template <typename TJournalString, typename TVariantMap, typename TVarKey, typename TPosition>
void _journalNextVariant(TJournalString & jString,
                         TVariantMap const & variantMap,
                         TVarKey const & varKey,
                         TPosition refCoordinate)
{
    switch(deltaType(varKey))
    {
        case DeltaType::DELTA_TYPE_SNP:
            _journalSnp(jString, refCoordinate, deltaSnp(variantMap, deltaPosition(varKey)));
            break;
        case DeltaType::DELTA_TYPE_DEL:
            _journalDel(jString, refCoordinate, deltaDel(variantMap, deltaPosition(varKey)));
            break;
        case DeltaType::DELTA_TYPE_INS:
            _journalIns(jString, refCoordinate, deltaIns(variantMap, deltaPosition(varKey)));
            break;
            // TODO(rmaerker): Add Case for INDEL
    }
}

// ----------------------------------------------------------------------------
// Function _doJournalBlock()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverage, typename TMapSpec, typename TSpec, typename TContextSize, typename TParallelTag>
inline void
_doJournalBlock(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverage, TMapSpec>, TSpec> & jst,
                TContextSize contextSize,
                Tag<TParallelTag> parallelTag = Serial())
{
    typedef JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverage, TMapSpec>, TSpec > TJst;
    typedef typename VariantData<TJst>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;
    typedef typename MappedDelta<TDeltaMap>::Type TMappedDelta;

    typedef typename Iterator<TDeltaCoverage, Standard>::Type TCoverageIterator;
    typedef typename Value<TDeltaCoverage>::Type TBitVec;
    typedef typename Iterator<TBitVec, Standard>::Type TBitVecIter;

    typedef typename JournalData<TJst>::Type TJournalSet;
    typedef typename Iterator<TJournalSet, Standard>::Type TJournalSetIter;
    typedef typename Size<TJournalSet>::Type TSize;

    // Define the block limits.
    TSize blockBegin = jst._activeBlock * jst._blockSize;
    TSize blockEnd = _min(length(keys(variantData(jst))), (jst._activeBlock + 1) * jst._blockSize);

    // Auxiliary variables.
    TDeltaMap & variantMap = variantData(jst);
    TJournalSet & journalSet = journalData(jst);
    String<int> _lastVisitedNodes;
    if (getBlockSize(jst) != TJst::REQUIRE_FULL_JOURNAL)
        resize(_lastVisitedNodes, journalSize(jst), -1, Exact());

    // Check whether there is enough space.
    SEQAN_ASSERT_EQ(length(journalSet), journalSize(jst));

    // Use parallel processing.
    // TODO(rmaerker): Consider more general Master-Worker design for parallelization?
    Splitter<TJournalSetIter> jSetSplitter(begin(journalSet, Standard()), end(journalSet, Standard()), parallelTag);
    SEQAN_OMP_PRAGMA(parallel for if (journalSize(jst) > 1000));
    for (unsigned jobId = 0; jobId < length(jSetSplitter); ++jobId)
    {

//        printf("Thread: %i of %i\n", jobId, omp_get_num_threads());
        unsigned jobBegin = jSetSplitter[jobId] - begin(journalSet, Standard());
        unsigned jobEnd = jSetSplitter[jobId + 1] - begin(journalSet, Standard());

        // Pre-processing: Update VPs for last block.
        if (getBlockSize(jst) != TJst::REQUIRE_FULL_JOURNAL)
            for (unsigned i = jobBegin; i < jobEnd; ++i)
            {
                clear(journalSet[i]);  // Reinitialize the journal strings.
                jst._blockVPOffset[i] += jst._activeBlockVPOffset[i];
            }

//        printf("Thread %i: jobBegin %i - jobEnd %i\n", jobId, jobBegin, jobEnd);

        TMapIterator itMapBegin = begin(variantMap, Standard());
        TMapIterator itMap = itMapBegin + blockBegin;
        TMapIterator itMapEnd = begin(variantMap, Standard()) + blockEnd;
        TCoverageIterator it = begin(deltaCoverageStore(variantMap), Standard()) + blockBegin;

        for (; itMap != itMapEnd; ++it, ++itMap)
        {
            TMappedDelta varKey = mappedDelta(variantMap, itMap - itMapBegin);

            TBitVecIter itVecBegin = begin(*it, Standard());
            TBitVecIter itVec = itVecBegin + jobBegin;
            TBitVecIter itVecEnd = begin(*it, Standard()) + jobEnd;
            for (;itVec != itVecEnd; ++itVec)
            {
                SEQAN_ASSERT_NOT(empty(host(journalSet[itVec - itVecBegin])));

                if (!(*itVec))
                    continue;

                // Store last visited node for current journaled string.
                if (getBlockSize(jst) != TJst::REQUIRE_FULL_JOURNAL)
                    _lastVisitedNodes[itVec - itVecBegin] = itMap - itMapBegin;
                _journalNextVariant(journalSet[itVec - itVecBegin], variantMap, varKey, *itMap);
            }
        }

        // Post-processing: Store VPs for current block.
        if (getBlockSize(jst) != TJst::REQUIRE_FULL_JOURNAL)
        {
            // Buffer the next branch nodes depending on the context size.
            for (unsigned i = jobBegin; i < jobEnd; ++i)
            {
                jst._activeBlockVPOffset[i] = length(journalSet[i]) - length(host(journalSet));
                if (_lastVisitedNodes[i] == -1)  // skip empty journal strings.
                    continue;

                TMappedDelta varKey = mappedDelta(variantMap, _lastVisitedNodes[i]);
                unsigned offset = keys(variantMap)[_lastVisitedNodes[i]] + contextSize;
                if (deltaType(varKey) == DeltaType::DELTA_TYPE_DEL)
                    offset += deltaDel(variantMap, deltaPosition(varKey));  // Adds the size of the deletion.
                // TODO(rmaerker): Add INDEL!

                if (itMapEnd != end(variantMap, Standard()) && offset >= *itMapEnd)
                {
                    // TODO(rmaerker): Check performance!
                    TMapIterator tmpMapIt = itMapEnd;
                    int localDiff = 0;
                    while (tmpMapIt != end(variantMap, Standard()) && *tmpMapIt <= offset + localDiff)
                    {
                        if (mappedCoverage(variantMap, tmpMapIt - itMapBegin)[i])
                        {
                            TMappedDelta varKey = mappedDelta(variantMap, tmpMapIt - itMapBegin);
                            _journalNextVariant(journalSet[i], variantMap, varKey, *tmpMapIt);

                            if (deltaType(varKey) == DeltaType::DELTA_TYPE_DEL)
                                localDiff += deltaDel(variantMap, deltaPosition(varKey));
                            else if (deltaType(varKey) == DeltaType::DELTA_TYPE_INS)
                                localDiff = _max(0, localDiff - static_cast<int>(length(deltaIns(variantMap, deltaPosition(varKey)))));
                        }
                        ++tmpMapIt;
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename THost>
inline void
setHost(JournaledStringTree<TDeltaMap, TSpec> & stringTree,
        THost & newHost)
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
// Function getBlockOffset()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec, typename TPosition>
inline typename MakeSigned<typename Size<JournaledStringTree<TDeltaMap, TSpec> const>::Type>::Type &
getBlockOffset(JournaledStringTree<TDeltaMap, TSpec> const & stringTree,
               TPosition const & pos)
{
    return stringTree._blockVPOffset[pos];
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

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
inline bool
requiresFullJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> const & stringTree)
{
    typedef JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> TJst;

    return stringTree._blockSize == TJst::REQUIRE_FULL_JOURNAL;
}

// ----------------------------------------------------------------------------
// Function requireJournal()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Deprecated!
template <typename TDeltaStore, typename TDeltaCoverageStore, typename TParallelTag>
void requireJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault> & stringTree,
                    Tag<TParallelTag> tag)
{
    SEQAN_ASSERT_NOT(empty(host(stringTree)));  // The host must be set before.
    if (stringTree._emptyJournal)
        _doJournalBlock(stringTree, tag);
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
void requireJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree)
{
    requireJournal(stringTree, Serial());
}

// ----------------------------------------------------------------------------
// Function journalNextBlock()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSize, typename TParallelTag>
bool journalNextBlock(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault> & stringTree,
                      TSize contextSize,
                      Tag<TParallelTag> tag)
{
    typedef JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, StringTreeDefault> TJst;

    if (stringTree._activeBlock < stringTree._numBlocks)
    {
        if (stringTree._blockSize == TJst::REQUIRE_FULL_JOURNAL && stringTree._emptyJournal)
        {
            _doJournalBlock(stringTree, contextSize, tag);
            stringTree._emptyJournal = false;
        }
        else if (stringTree._blockSize != TJst::REQUIRE_FULL_JOURNAL)
            _doJournalBlock(stringTree, contextSize, tag);

        ++stringTree._activeBlock;
        return true;
    }
    stringTree._activeBlock = 0;  // Reset the active block when all blocks are loaded.
    return false;
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TSize>
bool journalNextBlock(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree,
                      TSize contextSize)
{
    return journalNextBlock(stringTree, contextSize, Serial());
}

// ----------------------------------------------------------------------------
// Function reinitJournal()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TParallelTag>
void reinitJournal(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & jst,
                   Tag<TParallelTag> tag = Serial())
{
    typedef JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> TJst;
    typedef typename JournalData<TJst>::Type TJournaledSet;
    typedef typename Iterator<TJournaledSet, Standard>::Type TJournaledSetIterator;

    Splitter<TJournaledSetIterator> jsSplitter(begin(journalData(jst), Standard()), end(journalData(jst), Standard()), tag);
    SEQAN_OMP_PRAGMA(parallel for if (journalSize(jst) > 500))
    for (unsigned jobId = 0; jobId < length(jsSplitter); ++jobId)
    {
        for (TJournaledSetIterator threadIter = jsSplitter[jobId]; threadIter != jsSplitter[jobId+1]; ++threadIter)
            clear(*threadIter);
    }
//    clear(journalData(stringTree).strings);
//    clear(journalData(stringTree).limits);
//    journalData(stringTree).limitsValid = true;
//    createJournal(stringTree, tag);
}

// ----------------------------------------------------------------------------
// Function setBlockSize()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TPosition>
void setBlockSize(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> & stringTree,
                  TPosition const & newBlockSize)
{
    SEQAN_ASSERT_NOT(empty(stringTree._variantData));  // The delta map needs to be set before.

    stringTree._blockSize = newBlockSize;
    stringTree._activeBlock = 0;
    stringTree._numBlocks = static_cast<unsigned>(std::ceil(static_cast<double>(length(keys(variantData(stringTree)))) /
                                                            static_cast<double>(newBlockSize)));
    resize(stringTree._blockVPOffset, journalSize(stringTree), 0, Exact());
}

// ----------------------------------------------------------------------------
// Function getBlockSize()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec>
unsigned getBlockSize(JournaledStringTree<DeltaMap<TDeltaStore, TDeltaCoverageStore>, TSpec> const & stringTree)
{
    return stringTree._blockSize;
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

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_JOURNALED_STRING_TREE_DEFAULT_H_
