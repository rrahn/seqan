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
// Implements a data structure to store deltas efficiently.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DELTA_STORE_H_
#define EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DELTA_STORE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Struct DeltaType
// ----------------------------------------------------------------------------

struct DeltaType
{
    typedef size_t TValue;

    static const TValue MASK_DELTA_TYPE;
    static const TValue MASK_DELTA_POSITION;

    static const TValue DELTA_TYPE_SNP;
    static const TValue DELTA_TYPE_DEL;
    static const TValue DELTA_TYPE_INS;
    static const TValue DELTA_TYPE_INDEL;
};

// We make: 00 = SNP
//          01 = DEL
//          10 = INS
//          11 = INS_DEL -> INS_SNP can be replaced by INS_DEL.
const size_t DeltaType::MASK_DELTA_TYPE = 3ull << (BitsPerValue<size_t>::VALUE - 2);
const size_t DeltaType::MASK_DELTA_POSITION = ~MASK_DELTA_TYPE;
const size_t DeltaType::DELTA_TYPE_SNP = 0ull;
const size_t DeltaType::DELTA_TYPE_DEL = 1ull << (BitsPerValue<size_t>::VALUE - 2);
const size_t DeltaType::DELTA_TYPE_INS = 2ull << (BitsPerValue<size_t>::VALUE - 2);
const size_t DeltaType::DELTA_TYPE_INDEL = 3ull << (BitsPerValue<size_t>::VALUE - 2);

// ----------------------------------------------------------------------------
// Class DeltaStore
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
class DeltaStore
{
public:
    typedef String<TAlphabet> TSnpData;
    typedef String<TSize>   TDelData;
    typedef String<String<TAlphabet> > TInsData;
    typedef Pair<TSize, String<TAlphabet> > TInsDelDataValue;
    typedef String<TInsDelDataValue> TInsDelData;
    typedef typename DeltaType::TValue TValue;

    // TODO(rmaerker): Elaborate on these ideas!
    // Idea a) Use ConcatStringSet for insertions. Use as global insertion buffer for all journal sequences.
    // Idea b) Instead of storing all SNPs we store only a A,C,G,T, N -> AlphabetSize
    // Idea c) Instead of insertion buffer, we append the inserted strings to the reference and only use original nodes.

    String<TValue> _varDataMap;

    TSnpData    _snpData;
    TDelData    _delData;
    TInsData    _insData;
    TInsDelData _indelData;

    DeltaStore()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct Value<DeltaStore<TSize, TAlphabet> >
{
    typedef typename DeltaType::TValue Type;
};

template <typename TSize, typename TAlphabet>
struct Value<DeltaStore<TSize, TAlphabet> const>
{
    typedef typename DeltaType::TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct Reference<DeltaStore<TSize, TAlphabet> >
{
    typedef typename DeltaType::TValue & Type;
};

template <typename TSize, typename TAlphabet>
struct Reference<DeltaStore<TSize, TAlphabet> const>
{
    typedef typename DeltaType::TValue const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_SNP]
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename DeltaType::TValue>
struct DeltaValue;

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_SNP>
{
    typedef TAlphabet Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_SNP>
{
    typedef TAlphabet const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_DEL]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>
{
    typedef TSize Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_DEL>
{
    typedef TSize const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_INS]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>
{
    typedef String<TAlphabet> Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INS>
{
    typedef String<TAlphabet> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                   [DELTA_TYPE_INDEL]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INDEL>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>::Type TIns_;
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>::Type TDel_;
    typedef Pair<TDel_, TIns_> Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INDEL>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>::Type TIns_;
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>::Type TDel_;
    typedef Pair<TDel_, TIns_> const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function addSnpDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline void
addSnpDelta(DeltaStore<TSize, TAlphabet> & variantData, TAlphabet const & delta)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._snpData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._snpData, delta);
    appendValue(variantData._varDataMap, ((length(variantData._snpData) - 1) | DeltaType::DELTA_TYPE_SNP));
}

// ----------------------------------------------------------------------------
// Function addInsDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TStringSpec>
inline void
addInsDelta(DeltaStore<TSize, TAlphabet> & variantData,
            String<TAlphabet, TStringSpec> const & delta)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._insData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._insData, delta);
    appendValue(variantData._varDataMap, ((length(variantData._insData) - 1) | DeltaType::DELTA_TYPE_INS));
}

// ----------------------------------------------------------------------------
// Function addDelDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TSize2>
inline void
addDelDelta(DeltaStore<TSize, TAlphabet> & variantData,
            TSize2 const & delta)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._delData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._delData, delta);
    appendValue(variantData._varDataMap, ((length(variantData._delData) - 1) | DeltaType::DELTA_TYPE_DEL));
}

// ----------------------------------------------------------------------------
// Function addIndelDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TStringSpec, typename TSize2>
inline void
addIndelDelta(DeltaStore<TSize, TAlphabet> & variantData,
              Pair<TSize2, String<TAlphabet, TStringSpec> > const & delta)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._indelData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._indelData, delta);
    appendValue(variantData._varDataMap, ((length(variantData._indelData) - 1) | DeltaType::DELTA_TYPE_INDEL));
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename Reference<DeltaStore<TSize, TAlphabet> >::Type
value(DeltaStore<TSize, TAlphabet> & deltaStore,
      TPosition const & pos)
{
    return value(deltaStore._varDataMap, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename Reference<DeltaStore<TSize, TAlphabet> const>::Type
value(DeltaStore<TSize, TAlphabet> const & deltaStore,
      TPosition const & pos)
{
    return value(deltaStore._varDataMap, pos);
}

// ----------------------------------------------------------------------------
// Function deltaType()
// ----------------------------------------------------------------------------

inline DeltaType::TValue
deltaType(DeltaType::TValue const& val)
{
    return val & DeltaType::MASK_DELTA_TYPE;
}

// ----------------------------------------------------------------------------
// Function deltaPosition()
// ----------------------------------------------------------------------------

inline DeltaType::TValue
deltaPosition(DeltaType::TValue const& val)
{
    return val & DeltaType::MASK_DELTA_POSITION;
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._snpData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._snpData, pos);
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._delData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._delData, pos);
}

// ----------------------------------------------------------------------------
// Function deltaIns()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._insData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._insData, pos);
}


// ----------------------------------------------------------------------------
// Function deltaIndel()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._indelData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._indelData, pos);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DELTA_STORE_H_
