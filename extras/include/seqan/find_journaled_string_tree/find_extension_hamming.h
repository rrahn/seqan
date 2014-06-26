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
// Implements hamming distance based comparison of two strings.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag HammingPrefix
// ----------------------------------------------------------------------------

struct HammingPrefix_;
typedef Tag<HammingPrefix_> HammingPrefix;

// ----------------------------------------------------------------------------
// Tag HammingExtension
// ----------------------------------------------------------------------------

template <typename TTag = HammingSimple>
struct HammingExtension;

// ----------------------------------------------------------------------------
// Class Pattern                                                [HammingPrefix]
// ----------------------------------------------------------------------------

// NOTE(rmaerker): HammingPrefix is the same pattern as HammingSimple.
// We just need it to define a different iteration process.
template <typename TNeedle>
class Pattern<TNeedle, HammingPrefix> : Pattern<TNeedle, HammingSimple>
{
    typedef Pattern<TNeedle, HammingSimple> TSuper;

    Pattern() : TSuper()
    {}

    template <typename TNeedle2>
    Pattern(const TNeedle2 &ndl, int k = -1) : TSuper(ndl, k)
    {}
};

// ----------------------------------------------------------------------------
// Struct StateHammingSimple_
// ----------------------------------------------------------------------------

template <typename TSize>
struct StateHammingSimple_
{
    TSize  errorCount;
};

// ----------------------------------------------------------------------------
// Struct StateHammingPrefix_
// ----------------------------------------------------------------------------

template <typename TPattern_>
struct StateHammingPrefix_ : public StateHammingSimple_<int>
{
    typedef StateHammingSimple_<int>                    TSuper;
    typedef typename Host<TPattern>::Type               THost;
    typedef typename Iterator<THost, Standard>::Type    TNeedleIterator;

    TNeedleIterator currentIt;

    StateHammingPrefix_() : TSuper()
    {}
};

// ----------------------------------------------------------------------------
// Class FinderExtensionPoint                                   [HammingSimple]
// ----------------------------------------------------------------------------

template <typename TPattern_>
class FinderExtensionPoint<TPattern_,  HammingExtension<HammingSimple> > : FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_>             TSuper;
    typedef typename GetState<FinderExtensionPoint>::Type   TState;
    typedef typename Host<TPattern_>::Type                  THost;
    typedef typename Iterator<THost, Standard>::Type        TIterator;

    TState    _state;
    TIterator _ndlEnd;

    template <typename TErrors>
    FinderExtensionPoint(TPattern_ & pattern, TErrors errors) : TSuper(pattern)
    {
        init(*this, errors);
    }

    template <typename TResult, typename THystkIterator>
    inline void operator()(TResult & res, THystkIterator hystkIt)
    {
        _state.errorCount = 0;
        TIterator ndlIt = begin(host(getPattern(*this)), Standard());
        while(_state.errorCount <= getPattern(*this).maxDistance && ndlIt != _ndlEnd)
            if (*ndlIt++ != getValue(hystkIt++))
                ++_state.errorCount;
        res.i1 = _state.errorCount <= getPattern(*this).maxDistance;
    }
};

// ----------------------------------------------------------------------------
// Class FinderExtensionPoint                                   [HammingPrefix]
// ----------------------------------------------------------------------------

template <typename TPattern_>
class FinderExtensionPoint<TPattern_,  HammingExtension<HammingPrefix> > : FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_>             TSuper;
    typedef typename GetState<FinderExtensionPoint>::Type   TState;
    typedef typename Host<TPattern_>::Type                  THost;
    typedef typename Iterator<THost, Standard>::Type        TIterator;

    TState    _state;
    TIterator _ndlEnd;

    template <typename TErrors>
    FinderExtensionPoint(TPattern_ & pattern, TErrors errors) : TSuper(pattern)
    {
        init(*this, errors);
    }

    template <typename TResult, typename THystkIterator>
    inline void operator()(TResult & res, THystkIterator const & hystkIt)
    {
        if (*_state.currentIt != getValue(hystkIt))
            ++_state.errorCount;
        // Check if end of segment is reached.
        if (++_state.currentIt == _ndlEnd)
            res.i1 = _state.errorCount <= getPattern(*this).maxDistance;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition                               [Hamming]
// ----------------------------------------------------------------------------

// HammingSimple.
template <typename T>
struct ContextIteratorPosition<FinderExtensionPoint<T, HammingExtension<HammingSimple> > >
{
    typedef ContextPositionLeft Type;
};

// HammingPrefix.
template <typename T>
struct ContextIteratorPosition<FinderExtensionPoint<T, HammingExtension<HammingPrefix> > >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                    [Hamming]
// ----------------------------------------------------------------------------

// HammingPrefix.
template <typename T>
struct RequireFullContext<FinderExtensionPoint<T, HammingExtension<HammingPrefix> > > : False{};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                              [Hamming]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, HammingSimple> >
{
    typedef FinderExtensionPoint<Pattern<TNeedle, HammingSimple>, HammingExtension<HammingSimple> > Type;
};

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, HammingPrefix> >
{
    typedef FinderExtensionPoint<Pattern<TNeedle, HammingPrefix>, HammingExtension<HammingPrefix> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState                                              [Hamming]
// ----------------------------------------------------------------------------

template <typename T>
struct GetState<FinderExtensionPoint<T, HammingExtension<HammingSimple> > >
{
    typedef StateHammingSimple_<int> Type;
};

template <typename T>
struct GetState<FinderExtensionPoint<T const, HammingExtension<HammingSimple> > >
{
    typedef StateHammingSimple_<int> const Type;
};

template <typename T>
struct GetState<FinderExtensionPoint<T, HammingExtension<HammingPrefix> > >
{
    typedef StateHammingPrefix_<T> Type;
};

template <typename T>
struct GetState<FinderExtensionPoint<T const, HammingExtension<HammingPrefix> > >
{
    typedef StateHammingPrefix_<T> const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                                  [Hamming]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
inline typename GetState<FinderExtensionPoint<TPattern_, HammingExtension<TSpec> > >::Type &
getState(FinderExtensionPoint<TPattern_, HammingExtension<TSpec> > & extensionFunctor)
{
    return extensionFunctor._state;
}

template <typename TPattern_, typename TSpec>
inline typename GetState<FinderExtensionPoint<TPattern_, HammingExtension<TSpec> > const>::Type &
getState(FinderExtensionPoint<TPattern_, HammingExtension<TSpec> > const & extensionFunctor)
{
    return extensionFunctor._state;
}

// ----------------------------------------------------------------------------
// Function setState()                                                [Hamming]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec>
inline void
setState(FinderExtensionPoint<TPattern_, HammingExtension<TSpec> >  & extensionFunctor,
         typename GetState<FinderExtensionPoint<TPattern_, HammingExtension<TSpec> > >::Type const & state)
{
    extensionFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function initState()                                               [Hamming]
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline void
initState(FinderExtensionPoint<TPattern_, HammingExtension<HammingSimple> > & extensionFunctor)
{
    extensionFunctor._state.errorCount = 0;
}

template <typename TPattern_>
inline void
initState(FinderExtensionPoint<TPattern_, HammingExtension<HammingPrefix> > & extensionFunctor)
{
    extensionFunctor._state.errorCount = 0;
    extensionFunctor._state.currentIt = begin(host(getPattern(extensionFunctor)), Standard());
}

// ----------------------------------------------------------------------------
// Function init()                                                    [Hamming]
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TSpec, typename TMinScore>
inline void
init(FinderExtensionPoint<TPattern_, HammingExtension<TSpec> > & extension,
     TMinScore minScore)
{
    if (isInit(extension))
        return;

    setScoreLimit(getPattern(extension), minScore);
    initState(extension);
    extension._ndlEnd = end(host(getPattern(extension)), Standard());
    setInit(extension);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_
