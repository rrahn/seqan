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
// Implements the coverage store to hold the data of the coverage.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_COVERAGE_STORE_H_
#define EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_COVERAGE_STORE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DeltaCoverageStore
// ----------------------------------------------------------------------------

class DeltaCoverageStore
{
public:
    typedef String<bool, Packed<> > TBitVector;
    typedef Size<TBitVector>::Type TSize;

    TSize              _coverageSize;
    String<TBitVector> _coverageData;

    DeltaCoverageStore() : _coverageSize(0)
    {}

    template <typename TSize>
    DeltaCoverageStore(TSize const & newSize) : _coverageSize(newSize)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<DeltaCoverageStore>
{
    typedef String<bool, Packed<> > Type;
};

template <>
struct Value<DeltaCoverageStore const>
{
    typedef String<bool, Packed<> > const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <>
struct Reference<DeltaCoverageStore>
{
    typedef String<bool, Packed<> > & Type;
};

template <>
struct Reference<DeltaCoverageStore const>
{
    typedef String<bool, Packed<> > const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<DeltaCoverageStore>
{
    typedef typename Value<DeltaCoverageStore>::Type TValue_;
    typedef typename Size<TValue_>::Type Type;
};

template <>
struct Size<DeltaCoverageStore const>
{
    typedef typename Value<DeltaCoverageStore>::Type TValue_;
    typedef typename Size<TValue_>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<DeltaCoverageStore> :
    Size<DeltaCoverageStore>{};

template <>
struct Position<DeltaCoverageStore const> :
    Size<DeltaCoverageStore const>{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <>
struct Iterator<DeltaCoverageStore>
{
    typedef typename Value<DeltaCoverageStore>::Type TBitVector;
    typedef String<TBitVector> TCoverage;
    typedef typename Iterator<TCoverage, Standard>::Type Type;
};

template <>
struct Iterator<DeltaCoverageStore const>
{
    typedef typename Value<DeltaCoverageStore>::Type TBitVector;
    typedef String<TBitVector> TCoverage;
    typedef typename Iterator<TCoverage const, Standard>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

inline Iterator<DeltaCoverageStore, Standard>::Type
begin(DeltaCoverageStore & store, Standard const & /*tag*/)
{
    return begin(store._coverageData, Standard());
}

inline Iterator<DeltaCoverageStore const, Standard>::Type
begin(DeltaCoverageStore const & store, Standard const & /*tag*/)
{
    return begin(store._coverageData, Standard());
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

inline Iterator<DeltaCoverageStore, Standard>::Type
end(DeltaCoverageStore & store, Standard const & /*tag*/)
{
    return end(store._coverageData, Standard());
}

inline Iterator<DeltaCoverageStore const, Standard>::Type
end(DeltaCoverageStore const & store, Standard const & /*tag*/)
{
    return end(store._coverageData, Standard());
}

// ----------------------------------------------------------------------------
// Function addCoverage()
// ----------------------------------------------------------------------------

inline void
addCoverage(DeltaCoverageStore & store,
            Value<DeltaCoverageStore>::Type const & coverage)
{
    appendValue(store._coverageData, coverage);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TPosition>
inline typename Reference<DeltaCoverageStore>::Type
value(DeltaCoverageStore & store, TPosition const & pos)
{
    return value(store._coverageData, pos);
}

template <typename TPosition>
inline typename Reference<DeltaCoverageStore const>::Type
value(DeltaCoverageStore const & store, TPosition const & pos)
{
    return value(store._coverageData, pos);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TSize, typename TExpand>
inline void
resize(DeltaCoverageStore & store, TSize const & newSize, Tag<TExpand> const & tag)
{
    resize(store._coverageData, newSize, tag);
}

template <typename TSize>
inline void
resize(DeltaCoverageStore const & store, TSize const & newSize)
{
    typedef typename Value<DeltaCoverageStore>::Type TValue;
    typedef String<TValue> TData;

    resize(store, newSize, typename DefaultOverflowExplicit<TData>::Type());
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

inline Size<DeltaCoverageStore>::Type
length(DeltaCoverageStore const & store)
{
    return length(store._coverageData);
}

// ----------------------------------------------------------------------------
// Function setCoverageSize()
// ----------------------------------------------------------------------------

template <typename TSize>
inline void
setCoverageSize(DeltaCoverageStore & store, TSize newSize)
{
    store._coverageSize = newSize;
}

// ----------------------------------------------------------------------------
// Function coverageSize()
// ----------------------------------------------------------------------------

inline Size<DeltaCoverageStore>::Type
coverageSize(DeltaCoverageStore const & store)
{
    return store._coverageSize;
}

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_DELTA_COVERAGE_STORE_H_
