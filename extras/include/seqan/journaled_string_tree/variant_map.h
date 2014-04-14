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

template <typename TMap>
struct Keys;

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec = Default>
class DeltaMap
{
public:
    typedef TDeltaStore TDeltaStore_;
    typedef typename Keys<DeltaMap>::Type TKeys;

    TKeys               _keys;  // Key string containing the positions within the reference.
    TDeltaStore_         _deltaStore;
    TDeltaCoverageStore _deltaCoverageStore;

    DeltaMap() : _deltaStore(), _deltaCoverageStore()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Keys
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Keys<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >
{
    typedef typename Value<TDeltaStore>::Type TValue_;
    typedef String<TValue_> Type;
};

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Keys<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const >
{
    typedef typename Value<TDeltaStore>::Type TValue_;
    typedef String<TValue_> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Value<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >
{
    typedef typename Size<TDeltaStore>::Type Type;
};

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Value<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const >
{
    typedef typename Size<TDeltaStore>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Size<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >
{
    typedef typename Size<TDeltaStore>::Type Type;
};

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Size<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const >
{
    typedef typename Size<TDeltaStore>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Position<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >  :
    Size<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >{};

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct Position<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const >  :
    Size<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const>{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore>
struct Iterator<DeltaMap<TDeltaStore, TDeltaCoverageSore>, Standard>
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageSore> TDeltaMap_;
    typedef typename Keys<TDeltaMap_>::Type TKeys_;
    typedef typename Iterator<TKeys_, Standard>::Type Type;
};

template <typename TDeltaStore, typename TDeltaCoverageSore>
struct Iterator<DeltaMap<TDeltaStore, TDeltaCoverageSore> const, Standard>
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageSore> TDeltaMap_;
    typedef typename Keys<TDeltaMap_>::Type TKeys_;
    typedef typename Iterator<TKeys_ const, Standard>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultGetIteratorSpec
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore>
struct DefaultGetIteratorSpec<DeltaMap<TDeltaStore, TDeltaCoverageSore> const>
{
    typedef Standard Type;
};

// ----------------------------------------------------------------------------
// Metafunction MappedDelta
// ----------------------------------------------------------------------------

template <typename TMap>
struct MappedDelta;

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct MappedDelta<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >
{
    typedef typename Value<TDeltaStore>::Type Type;
};

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct MappedDelta<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const >
{
    typedef typename Value<TDeltaStore>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction MappedCoverage
// ----------------------------------------------------------------------------

template <typename TMap>
struct MappedCoverage;

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct MappedCoverage<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >
{
    typedef typename Value<TDeltaCoverageSore>::Type Type;
};

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
struct MappedCoverage<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const >
{
    typedef typename Value<TDeltaCoverageSore>::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageSore>, Standard>::Type
begin(DeltaMap<TDeltaStore, TDeltaCoverageSore> & variantStore,
      Standard const & /*tag*/)
{
    return begin(variantStore._keys, Standard());
}

template <typename TDeltaStore, typename TDeltaCoverageSore>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageSore> const, Standard>::Type
begin(DeltaMap<TDeltaStore, TDeltaCoverageSore> const & variantStore,
      Standard const & /*tag*/)
{
    return begin(variantStore._keys, Standard());
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageSore>, Standard>::Type
end(DeltaMap<TDeltaStore, TDeltaCoverageSore> & variantStore,
    Standard const & /*tag*/)
{
    return end(variantStore._keys, Standard());
}

template <typename TDeltaStore, typename TDeltaCoverageSore>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageSore> const, Standard>::Type
end(DeltaMap<TDeltaStore, TDeltaCoverageSore> const & variantStore,
    Standard const & /*tag*/)
{
    return end(variantStore._keys, Standard());
}

// ----------------------------------------------------------------------------
// Function deltaStore()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
inline TDeltaStore &
deltaStore(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore)
{
    return variantStore._deltaStore;
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
inline TDeltaStore const &
deltaStore(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore)
{
    return variantStore._deltaStore;
}

// ----------------------------------------------------------------------------
// Function deltaCoverageStore()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
inline DeltaCoverageStore &
deltaCoverageStore(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore)
{
    return variantStore._deltaCoverageStore;
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
inline DeltaCoverageStore const &
deltaCoverageStore(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore)
{
    return variantStore._deltaCoverageStore;
}

// ----------------------------------------------------------------------------
// Function keys()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
inline typename Keys<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >::Type &
keys(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore)
{
    return variantStore._keys;
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec>
inline typename Keys<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const>::Type &
keys(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore)
{
    return variantStore._keys;
}

// ----------------------------------------------------------------------------
// Function insert()                                                      [DEL]
// ----------------------------------------------------------------------------
// TODO(rmaerker): At the moment we just append the new value which is not clear by the name of this function. We need to adapt the function to do a binary search and then insert at the suggested position. For now the values must be presorted.

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec>, Standard>::Type
insert(DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec> & deltaMap,
       TKey const & key,
       typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_DEL>::Type const & delta)
{
    appendValue(keys(deltaMap), key);  // TODO(rmaerker): Adapt to binary search later.

    addDelDelta(deltaStore(deltaMap), delta);
    return end(deltaMap, Standard()) -1;
}

// ----------------------------------------------------------------------------
// Function insert()                                                      [SNP]
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec>, Standard>::Type
insert(DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec> & deltaMap,
       TKey const & key,
       typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_SNP>::Type const & delta)
{
    appendValue(keys(deltaMap), key);  // TODO(rmaerker): Adapt to binary search later.

    addSnpDelta(deltaStore(deltaMap), delta);
    return end(deltaMap, Standard()) -1;
}

// ----------------------------------------------------------------------------
// Function insert()                                                      [INS]
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec>, Standard>::Type
insert(DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec> & deltaMap,
       TKey const & key,
       typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_INS>::Type const & delta)
{
    appendValue(keys(deltaMap), key);  // TODO(rmaerker): Adapt to binary search later.

    addInsDelta(deltaStore(deltaMap), delta);
    return end(deltaMap, Standard()) -1;
}

// ----------------------------------------------------------------------------
// Function insert()                                                    [INDEL]
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec>, Standard>::Type
insert(DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec> & deltaMap,
       TKey const & key,
       typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_INDEL>::Type const & delta)
{
    appendValue(keys(deltaMap), key);  // TODO(rmaerker): Adapt to binary search later.

    addIndelDelta(deltaStore(deltaMap), delta);
    return end(deltaMap, Standard()) -1;
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TSpec, typename TKey,
          typename TDelta>
inline typename Iterator<DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec>, Standard>::Type
insert(DeltaMap<TDeltaStore, TDeltaCoverageStore, TSpec> & deltaMap,
       TKey const & key,
       TDelta const & delta,
       typename Value<TDeltaCoverageStore>::Type const & deltaCoverage)
{
    addCoverage(deltaCoverageStore(deltaMap), deltaCoverage);
    return insert(deltaMap, key, delta);
}

// ----------------------------------------------------------------------------
// Function mappedDelta()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPosition>
inline typename MappedDelta<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >::Type &
mappedDelta(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & deltaMap,
            TPosition const & pos)
{
   return value(deltaStore(deltaMap), pos);
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPosition>
inline typename MappedDelta<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const>::Type &
mappedDelta(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & deltaMap,
            TPosition const & pos)
{
   return value(deltaStore(deltaMap), pos);
}

// ----------------------------------------------------------------------------
// Function mappedCoverage()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPosition>
inline typename MappedCoverage<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> >::Type &
mappedCoverage(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore,
               TPosition const & pos)
{
    return value(deltaCoverageStore(variantStore), pos);
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPosition>
inline typename MappedCoverage<DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const>::Type &
mappedCoverage(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore,
               TPosition const & pos)
{
    return value(deltaCoverageStore(variantStore), pos);
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore,
         TPos const & pos)
{
    return deltaSnp(deltaStore(variantStore), pos);
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore,
         TPos const & pos)
{
    return deltaSnp(deltaStore(variantStore), pos);
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore,
         TPos const & pos)
{
    return deltaDel(deltaStore(variantStore), pos);
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore,
         TPos const & pos)
{
    return deltaDel(deltaStore(variantStore), pos);
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore,
         TPos const & pos)
{
    return deltaIns(deltaStore(variantStore), pos);
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore,
         TPos const & pos)
{
    return deltaIns(deltaStore(variantStore), pos);
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> & variantStore,
           TPos const & pos)
{
    return deltaIndel(deltaStore(variantStore), pos);
}

template <typename TDeltaStore, typename TDeltaCoverageSore, typename TSpec, typename TPos>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(DeltaMap<TDeltaStore, TDeltaCoverageSore, TSpec> const & variantStore,
           TPos const & pos)
{
    return deltaIndel(deltaStore(variantStore), pos);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_MAP_H_
