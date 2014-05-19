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
// Implements the iterator for a journaled index fibre.
// ==========================================================================
#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_JOURNALED_FIBRE_ITERATOR_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_JOURNALED_FIBRE_ITERATOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TJournaledFibre>
class Iter<TJournaledFibre, JournaledDeltaIndexFibreIterSpec>
{
public:
    typedef Iter<TJournaledFibre, JournaledDeltaIndexFibreIterSpec> TIterator;
    typedef typename GetJournalString_<TJournaledFibre>::Type TJournalString;
    typedef typename GetPartialSumManager_<TJournaledFibre>::Type TPartialSumManager;
    typedef typename Iterator<TJournalString, Standard>::Type TJournaledIter;

    TJournaledIter _hostIter;

    // TODO(rmaerker): Better to store the pointer to the underlying container?
    TPartialSumManager * _psMgrPtr;

    Iter() : _hostIter(), _psMgrPtr((TPartialSumManager*) 0)
    {}

    Iter(Iter const & other) : _hostIter(other._hostIter), _psMgrPtr(other._psMgrPtr)
    {}

    Iter(typename IterComplementConst<TIterator>::Type const & other) : _hostIter(other._hostIter), _psMgrPtr(other._psMgrPtr)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetPartialSumManager_
// ----------------------------------------------------------------------------

template <typename TJournaledIndexFibre_>
struct GetPartialSumManager_<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> >
{
    typedef typename GetPartialSumManager_<TJournaledIndexFibre_>::Type Type;
};

template <typename TJournaledIndexFibre_>
struct GetPartialSumManager_<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> const>
{
    typedef typename GetPartialSumManager_<TJournaledIndexFibre_ const>::Type Type;
};

// For String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

// TODO(rmaerker): fix doku!
///.Metafunction.Iterator.param.T:Spec.Journal String

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TTag>
struct Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, TTag>
{
    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> TJournaledIndexFibre_;
    typedef Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> Type;
};

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TTag>
struct Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, TTag>
{
    typedef JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const TJournaledIndexFibre_;
    typedef Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

// For Iter<TJournaledString, TJournaledStringIterSpec>
template <typename TJournaledIndexFibre_>
struct GetValue<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> >
{
    typedef typename GetValue<TJournaledIndexFibre_>::Type Type;
};

template <typename TJournaledIndexFibre_>
struct GetValue<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> const>
        : GetValue<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TJournaledIndexFibre_>
struct Value<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> >
{
    typedef typename Value<TJournaledIndexFibre_>::Type Type;
};

template <typename TJournaledIndexFibre_>
struct Value<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> const>
{
    typedef typename Value<TJournaledIndexFibre_ const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TJournaledIndexFibre_>
struct Reference<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> >
{
    typedef typename Reference<TJournaledIndexFibre_>::Type Type;
};

template <typename TJournaledIndexFibre_>
struct Reference<Iter<TJournaledIndexFibre_, JournaledDeltaIndexFibreIterSpec> const>
{
    typedef typename Reference<TJournaledIndexFibre_>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline
typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type
begin(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journalIndexFibre, Standard const &)
{
    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type TResult;
    TResult result;
    result._psMgrPtr = &value(journalIndexFibre._psMgrHolder);
    _initJournaledStringIterator(result._hostIter, journalIndexFibre._fibre);
    return result;
}

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline
typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, Standard>::Type
begin(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journalIndexFibre, Standard const &)
{
    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, Standard>::Type TResult;
    TResult result;
    result._psMgrPtr = &value(journalIndexFibre._psMgrHolder);
    _initJournaledStringIterator(result._hostIter, journalIndexFibre._fibre);
    return result;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline
typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type
end(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journalIndexFibre, Standard const &)
{
    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type TResult;
    TResult result;
    result._psMgrPtr = &value(journalIndexFibre._psMgrHolder);
    _initJournaledStringIteratorEnd(result._hostIter, journalIndexFibre._fibre);
    return result;
}

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager>
inline
typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, Standard>::Type
end(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journalIndexFibre, Standard const &)
{
    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, Standard>::Type TResult;
    TResult result;
    result._psMgrPtr = &value(journalIndexFibre._psMgrHolder);
    _initJournaledStringIteratorEnd(result._hostIter, journalIndexFibre._fibre);
    return result;
}

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPos>
inline
typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type
iter(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> & journalIndexFibre,
     TPos const & pos)
{
    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager>, Standard>::Type TResult;
    TResult result = begin(journalIndexFibre, Standard()) + pos;
    return result;
}

template <typename TValue, typename TJournaledSpec, typename TPartialSumManager, typename TPos>
inline
typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, Standard>::Type
iter(JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const & journalIndexFibre, TPos const & pos)
{
    typedef typename Iterator<JournalIndexFibre<TValue, TJournaledSpec, TPartialSumManager> const, Standard>::Type TResult;
    TResult result = begin(journalIndexFibre, Standard()) + pos;
    return result;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TJournalInfixFibre>
inline typename Reference<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> >::Type
value(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & me)
{
    typedef typename Reference<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> >::Type TReference;
    TReference ref(me);
    return ref;
}

template <typename TJournalInfixFibre>
inline typename Reference<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const>::Type
value(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & me)
{
    typedef typename Reference<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const>::Type TReference;
    TReference ref(me);
    return ref;
}

// assignValue

//template <typename TJournaledString, typename TJournalSpec, typename TValue>
//inline void
//assignValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & me,
//            TValue const & _value)
//{
//    SEQAN_CHECKPOINT;
//    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
//    assignValue(static_cast<TIterator const>(me), _value);
//}
//
//template <typename TJournaledString, typename TJournalSpec, typename TValue>
//inline void
//assignValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator,
//            TValue const & _value)
//{
//    SEQAN_CHECKPOINT;
//    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
//    typename Value<TIterator>::Type _temp_value = _value; //conversion
//    assignValue(*iterator._journalStringPtr, position(iterator), _temp_value);
//}

// moveValue

//template <typename TJournaledString, typename TJournalSpec, typename TValue>
//inline void
//moveValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & me,
//          TValue const & _value)
//{
//    SEQAN_CHECKPOINT;
//    // TODO(holtgrew): Copied from packed string. Actually, why no real move?
//    assignValue(me, _value);
//}
//
//template <typename TJournaledString, typename TJournalSpec, typename TValue>
//inline void
//moveValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & me,
//          TValue const & _value)
//{
//    SEQAN_CHECKPOINT;
//    // TODO(holtgrew): Copied from packed string. Actually, why no real move?
//    assignValue(me, _value);
//}

// valueConstruct

//template <typename TJournaledString, typename TJournalSpec>
//inline void
//valueConstruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & /*it*/)
//{
//    // TODO(holtgrew): Intentionally left blank? Alphabet elements must be default-constructable.
//}
//
//template <typename TJournaledString, typename TJournalSpec, typename TParam>
//inline void
//valueConstruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & it,
//               TParam const & param_)
//{
//    assignValue(it, param_);
//}
//
//template <typename TJournaledString, typename TJournalSpec, typename TParam>
//inline void
//valueConstruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & it,
//               TParam & param_,
//               Move const & /*tag*/)
//{
//    moveValue(it, param_);
//}
//
//// valueDestruct
//
//template <typename TJournaledString, typename TJournalSpec>
//inline void
//valueDestruct(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & /*it*/)
//{
//    // TODO(holtgrew): Intentionally left blank? Copied from packed string, leads to problems with non-POD contents!
//}

// position

template <typename TJournalInfixFibre>
inline typename Position<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const>::Type
position(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & iterator)
{
    return position(iterator._hostIter);
}

// ----------------------------------------------------------------------------
// Function _position()
// ----------------------------------------------------------------------------

// Returns the virtual position of the current iterator. Note, that the
// function position() should implement this behavior.
//template <typename TJournaledString, typename TJournalSpec>
//inline typename Position<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const>::Type
//_position(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
//{
//    if (atEnd(iterator._journalEntriesIterator))
//        return length(*iterator._journalStringPtr);
//
//    switch (value(iterator._journalEntriesIterator).segmentSource) {
//        case SOURCE_ORIGINAL:
//            return value(iterator._journalEntriesIterator).virtualPosition + iterator._currentHostIt - iterator._hostSegmentBegin;
//            break;
//        case SOURCE_PATCH:
//            return value(iterator._journalEntriesIterator).virtualPosition + iterator._currentInsertionBufferIt - iterator._insertionBufferSegmentBegin;
//            break;
//        default:
//            SEQAN_ASSERT_FAIL("Invalid segment source!");
//            return 0;
//    }
//}

// setPosition
//template <typename TJournaledString, typename TJournalSpec, typename TPosition>
//inline void
//setPosition(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & /*me*/,
//            TPosition /*pos_*/)
//{
//    SEQAN_CHECKPOINT;
//    // TODO(holtgrew): Implement me!
//    SEQAN_ASSERT_FAIL("Set position...");
//}

// getValue
template <typename TJournalInfixFibre>
inline typename GetValue<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const>::Type
getValue(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & me)
{
    if (me._hostIter._journalEntriesIterator->segmentSource == SOURCE_ORIGINAL)
    {
        return getValue(me._hostIter) + partialSum_(*me._psMgrPtr, getValue(me._hostIter));
    }
    else
        return getValue(me._hostIter);
}

template <typename TJournalInfixFibre>
inline typename GetValue<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> >::Type
getValue(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & me)
{
    if (me._hostIter._journalEntriesIterator->segmentSource == SOURCE_ORIGINAL)
    {
        return getValue(me._hostIter) + partialSum_(*me._psMgrPtr, getValue(me._hostIter));
    }
    else
        return getValue(me._hostIter);
}

template <typename TJournalInfixFibre>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> &
operator++(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator)
{
    ++iterator._hostIter;
    return iterator;
}

template <typename TJournalInfixFibre>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> &
operator++(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator, int /*postfix*/)
{
    Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> temp(iterator);
    ++iterator._hostIter;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator--()                                               [prefix]
// ----------------------------------------------------------------------------

template <typename TJournalInfixFibre>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> &
operator--(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator)
{
    --iterator._hostIter;
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator--()                                              [postfix]
// ----------------------------------------------------------------------------

template <typename TJournalInfixFibre>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> &
operator--(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator,
            int /*postfix*/)
{
    Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> temp(iterator);
    --iterator._hostIter;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

// TODO(rmaerker): We need a proxy class which is translated into value plus offset, if requested.
template <typename TJournalInfixFibre>
inline
typename Reference<TJournalInfixFibre>::Type
operator*(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator)
{
    return value(iterator);
}

template <typename TJournalInfixFibre>
inline
typename Reference<TJournalInfixFibre>::Type
operator*(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & iterator)
{
    return value(iterator);
}

// ----------------------------------------------------------------------------
// Funciton operator+=()
// ----------------------------------------------------------------------------

template <typename TJournalInfixFibre, typename TLen>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec>
operator+=(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator,
           TLen len_)
{
    iterator._hostIter += len_;
    return iterator;
}

template <typename TJournalInfixFibre>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec>
operator+(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator,
          typename Size<TJournalInfixFibre>::Type const & len)
{
    Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> temp(iterator);
    temp += len;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TJournalInfixFibre, typename TLen>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec>
operator-=(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator,
           TLen len_)
{
    iterator._hostIter -= len_;
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TJournalInfixFibre, typename TLen>
inline
Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec>
operator-(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> & iterator,
          typename Size<TJournalInfixFibre>::Type const & len)
{
    Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> temp(iterator);
    temp -= len;
    return temp;
}

template <typename TJournalInfixFibre>
inline
typename Difference<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> >::Type
operator-(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & it1,
          Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & it2)
{
    return it1._hostIter - it2._hostIter;
}

template <typename TJournalInfixFibre>
inline bool
operator==(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & a,
           Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & b)
{
    return a._hostIter == b._hostIter;
}

template <typename TJournalInfixFibre>
inline bool
operator==(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & a,
           typename IterComplementConst<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> >::Type const & b)
{
    return a._hostIter == b._hostIter;
}

template <typename TJournalInfixFibre>
inline bool
operator!=(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & a,
           Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & b)
{
    return !(a == b);
}

template <typename TJournalInfixFibre>
inline bool
operator!=(Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> const & a,
           typename IterComplementConst<Iter<TJournalInfixFibre, JournaledDeltaIndexFibreIterSpec> >::Type const & b)
{
    return !(a == b);
}

//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator<(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) < position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator<(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) < position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator<=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) <= position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator<=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) <= position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator>(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) > position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator>(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) > position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator>=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) >= position(b);
//}
//
//template <typename TJournaledString, typename TJournalSpec>
//inline
//bool
//operator>=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
//           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
//{
//    SEQAN_CHECKPOINT;
//    return position(a) >= position(b);
//}

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // SANDBOX_RMAERKER_INCLUDE_SEQAN_DELTA_INDEX_DELTA_INDEX_JOURNALED_FIBRE_ITERATOR_H_
