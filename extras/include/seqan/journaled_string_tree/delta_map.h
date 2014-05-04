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
// Implements the variant map. This is a facade that contains the delta
// store and the delta coverage store. It maps to each reference position
// the corresponding delta and coverage.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_MAP_H_
#define EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_MAP_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct SpecDeltaMapIterator_;

template <typename TMap>
struct GetMapValueString_;

template <typename TMap>
struct GetDeltaStore_;

template <typename TMap>
struct GetDeltaCoverageStore_;

template <typename TValue, typename TAlphabet, typename TSpec = Default>
class DeltaMap
{
public:

    typedef typename GetMapValueString_<DeltaMap>::Type TKeys;
    typedef typename GetDeltaStore_<DeltaMap>::Type TDeltaStore_;
    typedef typename GetDeltaCoverageStore_<DeltaMap>::Type TCoverageStore;

    TKeys            _keys;  // Key string containing the positions within the reference.
    TDeltaStore_     _deltaStore;
    TCoverageStore   _deltaCoverageStore;

    DeltaMap() : _keys(), _deltaStore(), _deltaCoverageStore()
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetMapValueString_
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetMapValueString_<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef String<TValue> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetMapValueString_<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef String<TValue> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetDeltaStore_
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaStore_<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename Size<TDeltaMap_>::Type TSize_;
    typedef DeltaStore<TSize_, TAlphabet> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaStore_<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename Size<TDeltaMap_>::Type TSize_;
    typedef DeltaStore<TSize_, TAlphabet> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetDeltaCoverageStore_
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaCoverageStore_<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef DeltaCoverageStore Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaCoverageStore_<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef DeltaCoverageStore const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct Value<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Value<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct Reference<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef TValue & Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Reference<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef TValue const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct Size<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename GetMapValueString_<TDeltaMap_>::Type TKeys_;

    typedef typename Size<TKeys_>::Type Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Size<DeltaMap<TValue, TAlphabet, TSpec> const > :
    Size<DeltaMap<TValue, TAlphabet, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct Position<DeltaMap<TValue, TAlphabet, TSpec> >  :
    Size<DeltaMap<TValue, TAlphabet, TSpec> >{};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Position<DeltaMap<TValue, TAlphabet, TSpec> const >  :
    Size<DeltaMap<TValue, TAlphabet, TSpec> const>{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, SpecDeltaMapIterator_> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_ const, SpecDeltaMapIterator_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultGetIteratorSpec
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct DefaultGetIteratorSpec<DeltaMap< TValue, TAlphabet, TSpec> >
{
    typedef SpecDeltaMapIterator_ Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct DefaultGetIteratorSpec<DeltaMap< TValue, TAlphabet, TSpec> const> :
    DefaultGetIteratorSpec<DeltaMap< TValue, TAlphabet, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename DeltaType::TValue TYPE>
struct DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, TYPE>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename TDeltaMap_::TDeltaStore_ TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_, TYPE>::Type Type;
};

template <typename TValue, typename TAlphabet, typename TSpec, typename DeltaType::TValue TYPE>
struct DeltaValue<DeltaMap<TValue, TAlphabet, TSpec> const, TYPE>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename TDeltaMap_::TDeltaStore_ TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_ const, TYPE>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaCoverage
// ----------------------------------------------------------------------------

template <typename TMap>
struct DeltaCoverage;

template <typename TValue, typename TAlphabet, typename TSpec>
struct DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef typename Value<DeltaCoverageStore>::Type Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef typename Value<DeltaCoverageStore>::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
find(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
     TKey const & key)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    SEQAN_ASSERT(!empty(deltaMap));

    TMapIterator it = std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), key);
    return (*it == key) ? it : end(deltaMap, Standard());
}

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
find(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
     TKey const & key)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> const TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    SEQAN_ASSERT(!empty(deltaMap));

    TMapIterator it = std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), key);
    return (*it == key) ? it : end(deltaMap, Standard());
}

// ----------------------------------------------------------------------------
// Function _insert()                                                     [DEL]
// ----------------------------------------------------------------------------
// TODO(rmaerker): At the moment we just append the new value which is not clear by the name of this function. We need to adapt the function to do a binary search and then insert at the suggested position. For now the values must be presorted.

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
       TKey const & key,
       TPosition const & insPos,
       typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_DEL>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);  // TODO(rmaerker): Adapt to binary search later.
    addDelDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function _insert()                                                     [SNP]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        TKey const & key,
        TPosition const & insPos,
        typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_SNP>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);  // TODO(rmaerker): Adapt to binary search later.
    addSnpDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function _insert()                                                     [INS]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        TKey const & key,
        TPosition const & insPos,
        typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_INS>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);
    addInsDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function _insert()                                                   [INDEL]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        TKey const & key,
        TPosition const & insPos,
        typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_INDEL>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);
    addIndelDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey,
          typename TDelta>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
       TKey const & key,
       TDelta const & delta,
       typename DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> >::Type const & deltaCoverage)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Position<TDeltaMap>::Type TPosition;

    TPosition insPos = std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), key) -
                       begin(deltaMap, Standard());
    addCoverage(deltaMap._deltaCoverageStore, deltaCoverage, insPos);
    _insert(deltaMap, key, insPos, delta);
    return iter(deltaMap, insPos, Standard());
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, Standard const & /*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type TIterator;
    TIterator tmp;
    _initBegin(tmp, deltaMap);
    return tmp;
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type TIterator;
    TIterator tmp;
    _initBegin(tmp, deltaMap);
    return tmp;
}

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
     TPos const & pos,
     Standard const & /*tag*/)
{
    return begin(deltaMap, Standard()) + pos;
}

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
     TPos const & pos,
     Standard const &/*tag*/)
{
    return begin(deltaMap, Standard()) + pos;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, Standard const &/*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type TIterator;
    TIterator tmp;
    _initEnd(tmp, deltaMap);
    return tmp;
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type TIterator;
    TIterator tmp;
    _initEnd(tmp, deltaMap);
    return tmp;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline void
clear(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap)
{
    clear(deltaMap._keys);
    clear(deltaMap._deltaStore);
    clear(deltaMap._deltaCoverageStore);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline bool
empty(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return empty(deltaMap._keys);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Size<DeltaMap<TValue, TAlphabet, TSpec> >::Type
length(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return length(deltaMap._keys);
}

// ----------------------------------------------------------------------------
// Function coverageSize()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Size<DeltaMap<TValue, TAlphabet, TSpec> >::Type
coverageSize(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return coverageSize(deltaMap._deltaCoverageStore);
}

// ----------------------------------------------------------------------------
// Function setCoverageSize()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TSize>
inline void
setCoverageSize(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, TSize const & size)
{
    return setCoverageSize(deltaMap._deltaCoverageStore, size);
}

// TODO(rmaerker): Implement erase
// TODO(rmaerker): Implement find
}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_MAP_H_
