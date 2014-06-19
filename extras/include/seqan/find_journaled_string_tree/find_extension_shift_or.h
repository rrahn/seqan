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
// Implements the shift or comparator.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_OR_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_OR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Functor Comparator()
// ----------------------------------------------------------------------------

template <typename TFinder>
class FinderExtensionPoint<TFinder, ShiftOr>
{
public:

    typedef typename GetPattern<TFinder>::Type TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    typedef unsigned int TWord;

    TWord       mask;
    TPattern    _state;
    bool        _isSmallNeedle;

    FinderExtensionPoint() : mask(), _state(), _isSmallNeedle(true)
    {}

    FinderExtensionPoint(TPattern & pattern)
    {
        init(*this, pattern);
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt const & haystackIt, BitAlgorithmSmallNeedle const &)
    {
        typedef typename Value<TNeedle>::Type TValue;

        _state.prefSufMatch[0] <<= 1;               //shift...
        _state.prefSufMatch[0] |= _state.bitMasks[ordValue(convert<TValue>(getValue(haystackIt)))];  //...or

        if (_state.prefSufMatch[0] & mask)  // If true, then there is no match.
            return;
        res.i1 = true;
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt const & haystackIt, BitAlgorithmLongNeedle const &)
    {
        typedef typename Value<TNeedle>::Type TValue;

        register TWord carry = 0;
        for(TWord block = 0; block < _state.blockCount; ++block)
        {
            bool newCarry = (_state.prefSufMatch[block] & (static_cast<TWord>(1) << (BitsPerValue<TWord>::VALUE - 1))) != 0;
            _state.prefSufMatch[block] <<= 1;
            _state.prefSufMatch[block] |= carry;
            carry = newCarry;
        }
        for(TWord block = 0; block < _state.blockCount; ++block)
            _state.prefSufMatch[block] |= _state.bitMasks[_state.blockCount * ordValue(convert<TValue>(getValue(haystackIt))) + block];
        if ((_state.prefSufMatch[_state.blockCount - 1] | mask) != static_cast<TWord>(~0))
        {
            res.i1 = true;
            return;
        }
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition                               [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct ContextIteratorPosition<FinderExtensionPoint<TFinder, ShiftOr> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                    [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct RequireFullContext<FinderExtensionPoint<TFinder, ShiftOr> > :
    False{};

// ----------------------------------------------------------------------------
// Metafunction GetState                                        [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct GetState<FinderExtensionPoint<TFinder, ShiftOr> >
{
    typedef FinderExtensionPoint<TFinder, ShiftOr> TExtensionfunctor;
    typedef typename TExtensionfunctor::TPattern Type;
};

template <typename TFinder>
struct GetState<FinderExtensionPoint<TFinder, ShiftOr> const>
{
    typedef FinderExtensionPoint<TFinder, ShiftOr> TExtensionfunctor;
    typedef typename TExtensionfunctor::TPattern const Type;
};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                              [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec>
struct RegisteredExtensionPoint<Finder_<TContainer, Pattern<TNeedle, ShiftOr>, Jst<TSpec> > >
{
    typedef Finder_<TContainer, Pattern<TNeedle, ShiftOr>, Jst<TSpec> > TFinder_;
    typedef FinderExtensionPoint<TFinder_, ShiftOr> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                                  [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
inline typename GetState<FinderExtensionPoint<TFinder, ShiftOr> >::Type &
getState(FinderExtensionPoint<TFinder, ShiftOr> & extensionFunctor)
{
    return extensionFunctor._state;
}

template <typename TFinder>
inline typename GetState<FinderExtensionPoint<TFinder, ShiftOr> const>::Type &
getState(FinderExtensionPoint<TFinder, ShiftOr> const & extensionFunctor)
{
    return extensionFunctor._state;
}

// ----------------------------------------------------------------------------
// Function initState()                                               [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
initState(FinderExtensionPoint<TFinder, ShiftOr> & extensionFunctor)
{
    _patternInit(extensionFunctor._state);
}

// ----------------------------------------------------------------------------
// Function setState()                                                [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
setState(FinderExtensionPoint<TFinder, ShiftOr>  & extensionFunctor,
                  typename GetState<FinderExtensionPoint<TFinder, ShiftOr> >::Type const & state)
{
    extensionFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function execute()
// ----------------------------------------------------------------------------

template <typename TResult, typename TFinder, typename TContextIter>
inline void
execute(TResult & res,
        FinderExtensionPoint<TFinder, ShiftOr> & extensionFunctor,
        TContextIter & contextIter)
{
    if (extensionFunctor._isSmallNeedle)
        extensionFunctor(res, contextIter, BitAlgorithmSmallNeedle());
    else
        extensionFunctor(res, contextIter, BitAlgorithmLongNeedle());
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
init(FinderExtensionPoint<TFinder, ShiftOr> & shiftOrFunctor,
     typename GetPattern<TFinder>::Type & pattern)
{
    typedef FinderExtensionPoint<TFinder, ShiftOr> TShiftOrFunctor;
    typedef typename TShiftOrFunctor::TWord TWord;

    _patternInit(pattern);
    shiftOrFunctor._state = pattern;

    if (length(host(pattern)) > BitsPerValue<TWord>::VALUE)
    {
        shiftOrFunctor._isSmallNeedle = false;
        shiftOrFunctor.mask = ~(static_cast<TWord>(1) << ((pattern.needleLength - 1) % BitsPerValue<TWord>::VALUE));
    }
    else
    {
        shiftOrFunctor._isSmallNeedle = true;
        shiftOrFunctor.mask = static_cast<TWord>(1) << (pattern.needleLength - 1);
    }
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_OR_H_
