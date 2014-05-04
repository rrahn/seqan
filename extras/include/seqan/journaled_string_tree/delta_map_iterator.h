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
// Implements the iterator over the delta map.
// ==========================================================================
#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TDeltaMap>
class Iter<TDeltaMap, SpecDeltaMapIterator_>
{
public:
    typedef typename GetMapValueString_<TDeltaMap>::Type TDeltaMapKeys;
    typedef typename Iterator<TDeltaMapKeys, Standard>::Type TMapIterator;

    typedef typename GetDeltaStore_<TDeltaMap>::Type TDeltaStore_;
    typedef typename Iterator<TDeltaStore_, Standard>::Type TDeltaStoreIter;

    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore_;
    typedef typename Iterator<TDeltaCoverageStore_, Standard>::Type TCoverageIter;

    TMapIterator    _mapIter;
    TDeltaStoreIter _deltaStoreIter;
    TCoverageIter   _deltaCoverageIter;

    Iter() : _mapIter(), _deltaStoreIter(), _deltaCoverageIter()
    {}

    // Copy Constructor
    Iter(Iter<TDeltaMap, SpecDeltaMapIterator_> const & other) : _mapIter(other._mapIter),
                                                                 _deltaStoreIter(other._deltaStoreIter),
                                                                 _deltaCoverageIter(other._deltaCoverageIter)
    {}

    Iter(typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & other) :
                    _mapIter(other._mapIter),
                    _deltaStoreIter(other._deltaStoreIter),
                   _deltaCoverageIter(other._deltaCoverageIter)
    {}

    // Assignment Operator
    Iter<TDeltaMap, SpecDeltaMapIterator_> &
    operator=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & other)
    {
        if (this != &other)
        {
            _mapIter = other._mapIter;
            _deltaStoreIter = other._deltaStoreIter;
            _deltaCoverageIter = other._deltaCoverageIter;
        }
        return *this;
    }

    Iter<TDeltaMap, SpecDeltaMapIterator_> &
    operator=(typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & other)
    {
        _mapIter = other._mapIter;
        _deltaStoreIter = other._deltaStoreIter;
        _deltaCoverageIter = other._deltaCoverageIter;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct Value<Iter<TDeltaMap, SpecDeltaMapIterator_> >
{
    typedef typename Value<TDeltaMap>::Type Type;
};

template <typename TDeltaMap>
struct Value<Iter<TDeltaMap, SpecDeltaMapIterator_> const >
{
    typedef typename Value<TDeltaMap const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct Reference<Iter<TDeltaMap, SpecDeltaMapIterator_> >
{
    typedef typename Reference<TDeltaMap>::Type Type;
};

template <typename TDeltaMap>
struct Reference<Iter<TDeltaMap, SpecDeltaMapIterator_> const >
{
    typedef typename Reference<TDeltaMap const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct Difference<Iter<TDeltaMap, SpecDeltaMapIterator_> >
{
    typedef Iter<TDeltaMap, SpecDeltaMapIterator_> TIter_;
    typedef typename TIter_::TMapIterator TMapIterator_;
    typedef typename Difference<TMapIterator_>::Type Type;
};

template <typename TDeltaMap>
struct Difference<Iter<TDeltaMap, SpecDeltaMapIterator_> const > :
    Difference<Iter<TDeltaMap, SpecDeltaMapIterator_> >{};

// ============================================================================
// Functions
// ============================================================================

template <typename TDeltaMap>
inline void
_initBegin(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter,
           TDeltaMap & map)
{
    iter._mapIter = begin(map._keys, Standard());
    iter._deltaStoreIter = begin(map._deltaStore, Standard());
    iter._deltaCoverageIter = begin(map._deltaCoverageStore, Standard());
}

//template <typename TDeltaMap>
//inline void
//_initBegin(Iter<TDeltaMap const, SpecDeltaMapIterator_> & iter,
//           TDeltaMap const & map)
//{
//    iter._mapIter = begin(map._keys, Standard());
//    iter._deltaStoreIter = begin(map._deltaStore, Standard());
//    iter._deltaCoverageIter = begin(map._deltaCoverageStore, Standard());
//}


template <typename TDeltaMap>
inline void
_initEnd(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter,
           TDeltaMap & map)
{
    iter._mapIter = end(map._keys, Standard());
    iter._deltaStoreIter = end(map._deltaStore, Standard());
    iter._deltaCoverageIter = end(map._deltaCoverageStore, Standard());
}

//template <typename TDeltaMap>
//inline void
//_initEnd(Iter<TDeltaMap const, SpecDeltaMapIterator_> & iter,
//           TDeltaMap const & map)
//{
//    iter._mapIter = end(map._keys, Standard());
//    iter._deltaStoreIter = end(map._deltaStore, Standard());
//    iter._deltaCoverageIter = end(map._deltaCoverageStore, Standard());
//}

template <typename TDeltaMap>
inline typename Value<TDeltaMap>::Type const &
value(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return value(iter._mapIter);
}

template <typename TDeltaMap>
inline typename Reference<TDeltaMap const>::Type
value(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return value(iter._mapIter);
}

// ----------------------------------------------------------------------------
// Function deltaType()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline DeltaType::TValue
deltaType(Iter<TDeltaMap, SpecDeltaMapIterator_> const& iter)
{
    return deltaType(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaPosition()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline DeltaType::TValue
deltaPosition(Iter<TDeltaMap, SpecDeltaMapIterator_> const& iter)
{
    return deltaPosition(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return deltaSnp(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return deltaSnp(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return deltaDel(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return deltaDel(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return deltaIns(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return deltaIns(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return deltaIndel(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return deltaIndel(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaCoverage()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename DeltaCoverage<TDeltaMap>::Type &
deltaCoverage(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return value(iter._deltaCoverageIter);
}

template <typename TDeltaMap>
inline typename DeltaCoverage<TDeltaMap>::Type &
deltaCoverage(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return value(iter._deltaCoverageIter);
}

// ----------------------------------------------------------------------------
// operator*
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Value<TDeltaMap>::Type const &
operator*(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    return value(iter);
}

template <typename TDeltaMap>
inline typename Value<TDeltaMap>::Type const &
operator*(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter)
{
    return value(iter);
}

// ----------------------------------------------------------------------------
// operator++
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline Iter<TDeltaMap, SpecDeltaMapIterator_> &
operator++(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    ++iter._mapIter;
    ++iter._deltaStoreIter;
    ++iter._deltaCoverageIter;
    return iter;
}

template <typename TDeltaMap>
inline Iter<TDeltaMap, SpecDeltaMapIterator_>
operator++(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter,  int /*postfix*/)
{
    Iter<TDeltaMap, SpecDeltaMapIterator_> temp(iter);
    ++iter;
    return temp;
}

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, SpecDeltaMapIterator_> &
operator+=(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter,  TSize const & len)
{
    iter._mapIter += len;
    iter._deltaStoreIter += len;
    iter._deltaCoverageIter += len;
    return iter;
}

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, SpecDeltaMapIterator_>
operator+(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter,  TSize const & len)
{
    Iter<TDeltaMap, SpecDeltaMapIterator_> temp(iter);
    temp += len;
    return temp;
}

template <typename TDeltaMap>
inline Iter<TDeltaMap, SpecDeltaMapIterator_> &
operator--(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter)
{
    --iter._mapIter;
    --iter._deltaStoreIter;
    --iter._deltaCoverageIter;
    return iter;
}

template <typename TDeltaMap>
inline Iter<TDeltaMap, SpecDeltaMapIterator_>
operator--(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter,  int /*postfix*/)
{
    Iter<TDeltaMap, SpecDeltaMapIterator_> temp(iter);
    --iter;
    return temp;
}

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, SpecDeltaMapIterator_> &
operator-=(Iter<TDeltaMap, SpecDeltaMapIterator_> & iter,  TSize const & len)
{
    iter._mapIter -= len;
    iter._deltaStoreIter -= len;
    iter._deltaCoverageIter -= len;
    return iter;
}

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, SpecDeltaMapIterator_>
operator-(Iter<TDeltaMap, SpecDeltaMapIterator_> const & iter,  TSize const & len)
{
    Iter<TDeltaMap, SpecDeltaMapIterator_> temp(iter);
    temp -= len;
    return temp;
}

template <typename TDeltaMap>
inline typename Difference<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type
operator-(Iter<TDeltaMap, SpecDeltaMapIterator_> const & lhs,
          Iter<TDeltaMap, SpecDeltaMapIterator_> const & rhs)
{
    return lhs._mapIter - rhs._mapIter;
}

template <typename TDeltaMap>
inline typename Difference<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type
operator-(Iter<TDeltaMap, SpecDeltaMapIterator_> const & lhs,
          typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & rhs)
{
    return lhs._mapIter - rhs._mapIter;
}

template <typename TDeltaMap>
inline bool
operator==(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           Iter<TDeltaMap, SpecDeltaMapIterator_> const & b)
{
    if (a._mapIter != b._mapIter)
        return false;
    return true;
}

template <typename TDeltaMap>
inline bool
operator==(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & b)
{
    typedef typename IterMakeConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type TConstIter;
    return static_cast<TConstIter>(a) == static_cast<TConstIter>(b);
}

template <typename TDeltaMap>
inline bool
operator!=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           Iter<TDeltaMap, SpecDeltaMapIterator_> const & b)
{
    return !(a == b);
}

template <typename TDeltaMap>
inline bool
operator!=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & b)
{
    return !(a == b);
}

template <typename TDeltaMap>
inline bool
operator<(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           Iter<TDeltaMap, SpecDeltaMapIterator_> const & b)
{
    return a._mapIter < b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator<(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
          typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & b)
{
    return a._mapIter < b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator<=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           Iter<TDeltaMap, SpecDeltaMapIterator_> const & b)
{
    return a._mapIter <= b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator<=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & b)
{
    return a._mapIter <= b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator>(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           Iter<TDeltaMap, SpecDeltaMapIterator_> const & b)
{
    return a._mapIter > b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator>(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & b)
{
    return a._mapIter > b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator>=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           Iter<TDeltaMap, SpecDeltaMapIterator_> const & b)
{
    return a._mapIter >= b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator>=(Iter<TDeltaMap, SpecDeltaMapIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaMap, SpecDeltaMapIterator_> >::Type const & b)
{
    return a._mapIter >= b._mapIter;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_
