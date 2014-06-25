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

template <typename TPattern_>
struct StateHammingSimple_
{
    typedef typename Host<TPattern_>::Type              THost;
    typedef typename Iterator<THost, Standard>::Type    TIterator;
    typedef typename Size<THost>::Type                  TSize;

    TSize       errorCount;
    TIterator   needleIt;
};

template <typename TPattern_>
class FinderExtensionPoint<TPattern_,  HammingSimple> : FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_>         TSuper;
    typedef typename Size<TPattern_>::Type              TSize;
//    typedef typename GetState<FinderExtensionPoint>::Type   TState;
    typedef typename Host<TPattern_>::Type              THost;
    typedef typename Iterator<THost, Standard>::Type    TIterator;

    TIterator _ndlBegin;
    TIterator _ndlEnd;

    template <typename TErrors>
    FinderExtensionPoint(TPattern_ & pattern, TErrors errors) : TSuper(pattern)
    {
        init(*this, errors);
    }

    template <typename TResult, typename THystkIterator>
    inline void operator()(TResult & res, THystkIterator hystkIt)
    {
        getPattern(*this).distance = 0;
        TIterator ndlIt = _ndlBegin;
        while(getPattern(*this).distance <= getPattern(*this).maxDistance && ndlIt != _ndlEnd)=
            if (*ndlIt++ != getValue(hystkIt++))
                ++getPattern(*this).distance;
        res.i1 = getPattern(*this).distance <= getPattern(*this).maxDistance;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition                               [Hamming]
// ----------------------------------------------------------------------------

template <typename T>
struct ContextIteratorPosition<FinderExtensionPoint<T, HammingSimple> >
{
    typedef ContextPositionLeft Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                    [Hamming]
// ----------------------------------------------------------------------------

//template <typename T>
//struct RequireFullContext<FinderExtensionPoint<T, HammingSimple> > : False{};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                              [Hamming]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, HammingSimple> >
{
    typedef FinderExtensionPoint<Pattern<TNeedle, HammingSimple>, HammingSimple> Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState                                              [Hamming]
// ----------------------------------------------------------------------------

//template <typename T>
//struct GetState<FinderExtensionPoint<T, HammingSimple> >
//{
//    typedef typename Size<T>::Type TSize_;
//    typedef StateHammingSimple_<TSize_> Type;
//};
//
//template <typename T>
//struct GetState<FinderExtensionPoint<T, HammingSimple> const>
//{
//    typedef typename Size<T>::Type TSize_;
//    typedef StateHammingSimple_<TSize_> const Type;
//};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                               [Pigeonhole]
// ----------------------------------------------------------------------------

//template <typename TPattern_, typename TSpec>
//inline typename GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > >::Type &
//getState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > & extensionFunctor)
//{
//    return extensionFunctor._state;
//}
//
//template <typename TPattern_, typename TSpec>
//inline typename GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > const>::Type &
//getState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > const & extensionFunctor)
//{
//    return extensionFunctor._state;
//}
//
//// ----------------------------------------------------------------------------
//// Function setState()                                             [Pigeonhole]
//// ----------------------------------------------------------------------------
//
//template <typename TPattern_, typename TSpec>
//inline void
//setState(FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> >  & extensionFunctor,
//         typename GetState<FinderExtensionPoint<TPattern_, Pigeonhole<TSpec> > >::Type const & state)
//{
//    extensionFunctor._state = state;
//}

// ----------------------------------------------------------------------------
// Function initState()
// ----------------------------------------------------------------------------

//template <typename TPattern_>
//inline void
//initState(FinderExtensionPoint<TPattern_, HammingSimple> & extensionFunctor)
//{
//    extensionFunctor._state.errorCount = 0;
//    extensionFunctor._state.needleIt = begin(host(getPattern(extensionFunctor)), Standard());
//}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TMinScore>
inline void
init(FinderExtensionPoint<TPattern_, HammingSimple> & extension,
     TMinScore minScore)
{
    if (isInit(extension))
        return;

    setScoreLimit(getPattern(extension), minScore);
    extension._ndlBegin = begin(host(getPattern(extension)), Standard());
    extension._ndlEnd = end(host(getPattern(extension)), Standard());
    setInit(extension);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_
