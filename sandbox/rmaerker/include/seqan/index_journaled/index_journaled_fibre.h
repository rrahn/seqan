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
// Implements the adapted version of the journaled string, which additionally
// stores the offset of values.
// ==========================================================================
#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_JOURNALED_FIBRE_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_JOURNALED_FIBRE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
class JournalIndexFibre
{
public:
    typedef String<TValue, TJournaledSpec> TFibre;
    typedef Holder<TPartialSumManager> TPSMgrHolder;
//    typedef TPartialSumManager* TPSMgrPtr;


    TFibre _fibre;
    TPSMgrHolder _psMgrHolder;

    // Default c'tor.
    JournalIndexFibre() : _fibre(), _psMgrHolder()
    {}

    // TODO(rmaerker): Why do we need two versions for const and non const journal string.
    // Standard c'tor with host string.
    template <typename TSpec>
    JournalIndexFibre(String<TValue, TSpec> & hostString) :
        _fibre(hostString),
        _psMgrHolder()
    {}

    template <typename TSpec>
    JournalIndexFibre(String<TValue, TSpec> const & hostString) :
        _fibre(hostString),
        _psMgrHolder()

    {}

    // Copy c'tor
    JournalIndexFibre(JournalIndexFibre & other) : _fibre(other._fibre), _psMgrHolder(other._psMgrHolder)
    {}

    JournalIndexFibre(JournalIndexFibre const & other) : _fibre(other._fibre), _psMgrHolder(other._psMgrHolder)
    {}

    // Assignment operator.
    JournalIndexFibre &
    operator=(JournalIndexFibre & other)
    {
        if (this != &other)
        {
            _fibre = other._fibre;
            _psMgrHolder = other._psMgrHolder;
        }
        return *this;
    }

    JournalIndexFibre &
    operator=(JournalIndexFibre const & other)
    {
        if (this != &other)
        {
            _fibre = other._fibre;
            _psMgrHolder = other._psMgrHolder;
        }
        return *this;
    }
};



// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetPartialSumManager_
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct GetPartialSumManager_<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >
{
    typedef TPartialSumManager Type;
};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct GetPartialSumManager_<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>
{
    typedef TPartialSumManager const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct Host<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >
{
    typedef String<TValue, TJournaledSpec> TJournalString;
    typedef typename Host<TJournalString>::Type Type;
};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct Host<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>
{
    typedef String<TValue, TJournaledSpec> const TJournalString;
    typedef typename Host<TJournalString>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct Value<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> > :
    Value<String<TValue, TJournaledSpec> > {};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct Value<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const> :
    Value<String<TValue, TJournaledSpec> const> {};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct GetValue<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> > :
    Value<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> > {};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct GetValue<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const> :
    Value<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const> {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct Reference<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >
{
    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> TJournalFibre_;
    typedef typename Iterator<TJournalFibre_, Standard>::Type TIter_;
    typedef Proxy<IteratorProxy<TIter_> > TProxy_;
    typedef TProxy_ Type;
};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct Reference<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>
{
    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const TJournalFibre_;
    typedef typename Iterator<TJournalFibre_, Standard>::Type TIter_;
    typedef Proxy<IteratorProxy<TIter_> > TProxy_;
    typedef TProxy_ Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetJournaledString_
// ----------------------------------------------------------------------------

template <typename TJournalIndexFibre_>
struct GetJournalString_;

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct GetJournalString_<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >
{
    typedef String<TValue, TJournaledSpec> Type;
};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
struct GetJournalString_<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>
{
    typedef String<TValue, TJournaledSpec> const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Consider special ProxyIterator which stores the ProxyIterator and the offset.
// TODO(rmaerker): Write docu.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPosition>
inline typename Reference<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >::Type
value(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre, TPosition const & pos)
{
    //    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type TIter;
    //    typedef typename Reference<TIter>::Type TRefIter;

//    static_cast<Nothing>(TIter());
//    TIter it = iter(journaledIndexFibre, pos);
//    static_cast<Nothing>(it);
//    static_cast<Nothing>(*it);
//    TRefIter ref = *it;
//    static_cast<Nothing>(ref);
//    TRefIter ref2 = *iter(journaledIndexFibre, pos);
    return *iter(journaledIndexFibre, pos);
}

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPosition>
inline typename Reference<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>::Type
value(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journaledIndexFibre, TPosition const & pos)
{
    return *iter(journaledIndexFibre, pos);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Write docu.
// Returns the actual value which corresponds to the pointed value and the offset to this position.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPosition>
inline typename GetValue<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>::Type
getValue(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journaledIndexFibre, TPosition const & pos)
{
    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const TIndexFibre;
    typedef typename TIndexFibre::TFibre TJournalString;
    typedef typename TJournalString::TJournalEntry TJournalEntry;
    typedef typename Position<TJournalString>::Type TPos;

    TJournalEntry entry = findJournalEntry(journaledIndexFibre._fibre._journalEntries, pos);
    TPos relativePos = pos - entry.virtualPosition;

    if (entry.segmentSource == SOURCE_ORIGINAL)
    {
        TValue val = getValue(value(journaledIndexFibre._fibre._holder), entry.physicalPosition + relativePos);
        return val + partialSum_(psManager(journaledIndexFibre), val);
    }
    else
    {
        return getValue(journaledIndexFibre._fibre._insertionBuffer, entry.physicalPosition + relativePos);
    }
//    return  val + partialSum_(*journaledIndexFibre._psMgrPtr, val);
}

//template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPosition>
//inline typename GetValue<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >::Type
//getValue(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre, TPosition const & pos)
//{
//    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const TIndexFibre;
//    typedef typename TIndexFibre::TFibre TJournalString;
//    typedef typename TJournalString::TJournalEntry TJournalEntry;
//    typedef typename Position<TJournalString>::Type TPos;
//
//    TJournalEntry entry = findJournalEntry(journaledIndexFibre._fibre._journalEntries, pos);
//    TPos relativePos = pos - entry.virtualPosition;
//
//    if (entry.segmentSource == SOURCE_ORIGINAL)
//    {
//        TValue val = getValue(value(journaledIndexFibre._fibre._holder), entry.physicalPosition + relativePos);
//        return val + partialSum_(psManager(journaledIndexFibre), val);
//    }
//    else
//    {
//        return getValue(journaledIndexFibre._fibre._insertionBuffer, entry.physicalPosition + relativePos);
//    }
////    return  val + partialSum_(*journaledIndexFibre._psMgrPtr, val);
//}

//template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPosition>
//inline typename GetValue<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >::Type
//getValue(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre, TPosition const & pos)
//{
//    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> TDeltaIndex;
//    typedef typename SAValue<TDeltaIndex>::Type TSAValue;
//
//    TSAValue val = getValue(journaledIndexFibre._fibre, pos);
//    return  val + partialSum_(*journaledIndexFibre._psMgrPtr, val);
//}

// ----------------------------------------------------------------------------
// Function createPSManager()
// ----------------------------------------------------------------------------

// Creates a copy of the psManager and stores its copy in the holder.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline void
createPSManager(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre,
                TPartialSumManager const & psMgr)
{
    create(journaledIndexFibre._psMgrHolder, psMgr);
}

// ----------------------------------------------------------------------------
// Function setPSManager()
// ----------------------------------------------------------------------------

// Stores a reference to the underlying psMgr.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline void
setPSManager(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre,
                TPartialSumManager & psMgr)
{
    setValue(journaledIndexFibre._psMgrHolder, psMgr);
}

// Stores a reference to the underlying psMgr.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline void
setPSManager(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre,
                TPartialSumManager const & psMgr)
{
    setValue(journaledIndexFibre._psMgrHolder, psMgr);
}

// ----------------------------------------------------------------------------
// Function psManager()
// ----------------------------------------------------------------------------

// Stores a reference to the underlying psMgr.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline TPartialSumManager &
psManager(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre)
{
    return value(journaledIndexFibre._psMgrHolder);
}

// Stores a reference to the underlying psMgr.
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline TPartialSumManager const &
psManager(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journaledIndexFibre)
{
    return value(journaledIndexFibre._psMgrHolder);
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Fix docu.
/**
.Function.setHost:
..class:Spec.Journaled String
..param.object.type:Spec.Journaled String
..param.host.type:Class.String
..include:seqan/sequence_journaled.h
*/
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TSequence2>
inline void
setHost(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre, TSequence2 & host)
{
    setValue(journaledIndexFibre._fibre._holder, host);
    journaledIndexFibre._fibre._length = length(host);
    reinit(journaledIndexFibre._fibre, length(host));
}

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Fix docu.
/**
.Function.host:
..class:Spec.Journaled String
..param.object.type:Spec.Journaled String
..include:seqan/sequence_journaled.h
*/
template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline
typename Host<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> >::Type &
host(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre)
{
    return value(journaledIndexFibre._fibre._holder);
}

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline
typename Host<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const>::Type &
host(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journaledIndexFibre)
{
    return value(journaledIndexFibre._fibre._holder);
}

// ----------------------------------------------------------------------------
// Function reinit()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline void
reinit(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre)
{
    reinit(journaledIndexFibre._fibre);
    clear(value(journaledIndexFibre._psMgrHolder));  // TODO(rmaerker): Other tables that rely on same manager need to be informed if its state is changed.
}

// ----------------------------------------------------------------------------
// Function lenght()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline typename Size<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const >::Type
length(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journaledIndexFibre)
{
    return length(journaledIndexFibre._fibre);
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

//template <typename TJournaledString, typename TPos>
//inline void
//insert(JournaledIndexTable_<TJournaledString, Simple> & journaledIndexTable,
//       TPos pos,
//       typename Host<TJournaledString>::Type const & seq)
//{
//    insert(journaledIndexTable._journaledTable, pos, seq);
//
//    journaledString._length += length(seq);
//    TPos beginPos = length(journaledString._insertionBuffer);
//    append(journaledString._insertionBuffer, seq);
//    recordInsertion(journaledString._journalEntries, pos, beginPos, length(seq));
//}

// ----------------------------------------------------------------------------
// Function insertValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TInsertPos, typename TValue2>
inline void
insertValue(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre,
            TInsertPos pos,
            TValue2 const & value)
{
    insertValue(journaledIndexFibre._fibre, pos, value);
}

// ----------------------------------------------------------------------------
// Function erase()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TBeginPos, typename TEndPos>
inline void
erase(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre,
      TBeginPos pos,
      TEndPos posEnd)
{
    erase(journaledIndexFibre._fibre, pos, posEnd);
}

// ----------------------------------------------------------------------------
// Function erase()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPos>
inline void
erase(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre,
      TPos pos)
{
    erase(journaledIndexFibre, pos, pos + 1);
}


template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline void
_printDebug(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journaledIndexFibre)
{
    std::cout << "##################### DEBUG OUTPUT #####################" << std::endl;
    std::cout << "The Host: " << host(journaledIndexFibre) << std::endl;
    std::cout << "The Journal Entries: " << journaledIndexFibre._fibre._journalEntries << std::endl;
    std::cout << "The Insertion Buffer: " << journaledIndexFibre._fibre._insertionBuffer << std::endl;
    std::cout << "The PartialSum Mgr: " << value(journaledIndexFibre._psMgrHolder)._partialSumTable << std::endl;
    std::cout << "######################### CLOSE #########################" << std::endl;
}

}  // namespace seqan

#endif  // SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_JOURNALED_FIBRE_H_
