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
// Implements the shift and search.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Functor FinderExtensionPoint
// ----------------------------------------------------------------------------

template <typename TFinder>
class FinderExtensionPoint<TFinder, ShiftAnd >
{
public:

    typedef typename GetPattern<TFinder>::Type TPattern;
    typedef typename Needle<TPattern>::Type TNeedle;
    typedef unsigned int TWord;

    TPattern  _state;
    TWord     _compare;
    bool      _isSmallNeedle;

    FinderExtensionPoint() : _state(), _compare(), _isSmallNeedle()
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

        _state.prefSufMatch[0] = ((_state.prefSufMatch[0] << 1) | static_cast<TWord>(1)) &
                                       _state.bitMasks[ordValue(convert<TValue>(getValue(haystackIt)))];
        if ((_state.prefSufMatch[0] & _compare) != 0)
            res.i1 = true;
    }

    template <typename TResult, typename THystkIt>
    inline void
    operator()(TResult & res, THystkIt const & haystackIt, BitAlgorithmLongNeedle const &)
    {
        typedef typename Value<TNeedle>::Type TValue;

        register TWord carry = 1; // Only set when going into the position.
        // Is set in each iteration.
        for(TWord block = 0; block < _state.blockCount; ++block)
        {
            bool newCarry = (_state.prefSufMatch[block] & (static_cast<TWord>(1) << (BitsPerValue<TWord>::VALUE - 1))) != 0;
            _state.prefSufMatch[block] <<= 1;
            _state.prefSufMatch[block] |= carry;
            carry = newCarry;
        }
        for(TWord block = 0; block < _state.blockCount; ++block)
            _state.prefSufMatch[block] &= _state.bitMasks[_state.blockCount * ordValue(convert<TValue>(getValue(haystackIt))) + block];
        if ((_state.prefSufMatch[_state.blockCount-1] & _compare) != 0)
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
// Metafunction ContextIteratorPosition                              [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct ContextIteratorPosition<FinderExtensionPoint<TFinder, ShiftAnd> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                   [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TFinder>
struct RequireFullContext<FinderExtensionPoint<TFinder, ShiftAnd> > :
    False{};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TFinder>
struct GetState<FinderExtensionPoint<TFinder, ShiftAnd> >
{
    typedef FinderExtensionPoint<TFinder, ShiftAnd> TExtensionfunctor;
    typedef typename TExtensionfunctor::TPattern Type;
};

template <typename TFinder>
struct GetState<FinderExtensionPoint<TFinder, ShiftAnd> const>
{
    typedef FinderExtensionPoint<TFinder, ShiftAnd> TExtensionfunctor;
    typedef typename TExtensionfunctor::TPattern const Type;
};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                             [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec>
struct RegisteredExtensionPoint<Finder_<TContainer, Pattern<TNeedle, ShiftAnd>, Jst<TSpec> > >
{
    typedef Finder_<TContainer, Pattern<TNeedle, ShiftAnd>, Jst<TSpec> > TFinder_;
    typedef FinderExtensionPoint<TFinder_, ShiftAnd> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                                 [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TFinder>
inline typename GetState<FinderExtensionPoint<TFinder, ShiftAnd> >::Type &
getState(FinderExtensionPoint<TFinder, ShiftAnd> & extensionFunctor)
{
    return extensionFunctor._state;
}

template <typename TFinder>
inline typename GetState<FinderExtensionPoint<TFinder, ShiftAnd> const>::Type &
getState(FinderExtensionPoint<TFinder, ShiftAnd> const & extensionFunctor)
{
    return extensionFunctor._state;
}

// ----------------------------------------------------------------------------
// Function setState()                                               [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
setState(FinderExtensionPoint<TFinder, ShiftAnd>  & extensionFunctor,
         typename GetState<FinderExtensionPoint<TFinder, ShiftAnd> >::Type const & state)
{
    extensionFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function initState()                                               [ShiftOr]
// ----------------------------------------------------------------------------

template <typename TFinder>
inline void
initState(FinderExtensionPoint<TFinder, ShiftAnd> & extensionFunctor)
{
    _patternInit(extensionFunctor._state);
}

// ----------------------------------------------------------------------------
// Function execute()
// ----------------------------------------------------------------------------

template <typename TResult, typename TFinder, typename TContextIter>
inline void
execute(TResult & res,
        FinderExtensionPoint<TFinder, ShiftAnd> & extensionFunctor,
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
init(FinderExtensionPoint<TFinder, ShiftAnd> & shiftAndFunctor,
     typename GetPattern<TFinder>::Type & pattern)
{
    typedef FinderExtensionPoint<TFinder, ShiftAnd> TShiftAndFunctor;
    typedef typename TShiftAndFunctor::TWord TWord;

    _patternInit(pattern);
    shiftAndFunctor._state = pattern;
    if (length(host(pattern)) > BitsPerValue<TWord>::VALUE)
    {
        shiftAndFunctor._isSmallNeedle = false;
        shiftAndFunctor._compare = static_cast<TWord>(1) << ((pattern.needleLength - 1) % BitsPerValue<TWord>::VALUE);
    }
    else
    {
        shiftAndFunctor._isSmallNeedle = true;
        shiftAndFunctor._compare = static_cast<TWord>(1) << (pattern.needleLength - 1);
    }
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_
