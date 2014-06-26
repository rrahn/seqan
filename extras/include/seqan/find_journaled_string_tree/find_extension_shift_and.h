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

template <typename TPattern_>
class FinderExtensionPoint<TPattern_, ShiftAnd> : public FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_> TSuper;
    typedef typename Host<TPattern_>::Type TNeedle;
    typedef typename GetState<FinderExtensionPoint>::Type TState;
    typedef unsigned int TWord;

    TState      _state;
    TWord       _compare;
    bool        _isSmallNeedle;

    FinderExtensionPoint(TPattern_ & pattern) : TSuper(pattern)
    {
        init(*this);
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

template <typename TPattern_>
struct ContextIteratorPosition<FinderExtensionPoint<TPattern_, ShiftAnd> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                   [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TPattern_>
struct RequireFullContext<FinderExtensionPoint<TPattern_, ShiftAnd> > :
    False{};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TPattern_>
struct GetState<FinderExtensionPoint<TPattern_, ShiftAnd> >
{
    typedef TPattern_ Type;
};

template <typename TPattern_>
struct GetState<FinderExtensionPoint<TPattern_, ShiftAnd> const>
{
    typedef TPattern_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                             [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, ShiftAnd> >
{
    typedef FinderExtensionPoint<Pattern<TNeedle, ShiftAnd>, ShiftAnd> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState                                        [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline typename GetState<FinderExtensionPoint<TPattern_, ShiftAnd> >::Type &
getState(FinderExtensionPoint<TPattern_, ShiftAnd> & extensionFunctor)
{
    return extensionFunctor._state;
}

template <typename TPattern_>
inline typename GetState<FinderExtensionPoint<TPattern_, ShiftAnd> const>::Type &
getState(FinderExtensionPoint<TPattern_, ShiftAnd> const & extensionFunctor)
{
    return extensionFunctor._state;
}

// ----------------------------------------------------------------------------
// Function setState()                                      [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline void
setState(FinderExtensionPoint<TPattern_, ShiftAnd>  & extensionFunctor,
         typename GetState<FinderExtensionPoint<TPattern_, ShiftAnd> >::Type const & state)
{
    extensionFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function initState()
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline void
initState(FinderExtensionPoint<TPattern_, ShiftAnd> & extension)
{
    extension._state = extension._pattern;
}

// ----------------------------------------------------------------------------
// Function execute()
// ----------------------------------------------------------------------------

template <typename TResult, typename TPattern_, typename TContextIter>
inline void
execute(TResult & res,
        FinderExtensionPoint<TPattern_, ShiftAnd> & extensionFunctor,
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

template <typename TPattern_>
inline void
init(FinderExtensionPoint<TPattern_, ShiftAnd> & extension)
{
    typedef FinderExtensionPoint<TPattern_, ShiftAnd> TShiftAndFunctor;
    typedef typename TShiftAndFunctor::TWord TWord;

    if (isInit(extension))
        return;

    _patternInit(getPattern(extension));
    initState(extension);

    if (contextSize(extension) > BitsPerValue<TWord>::VALUE)
    {
        extension._isSmallNeedle = false;
        extension._compare = static_cast<TWord>(1) << ((getPattern(extension).needleLength - 1) % BitsPerValue<TWord>::VALUE);
    }
    else
    {
        extension._isSmallNeedle = true;
        extension._compare = static_cast<TWord>(1) << (getPattern(extension).needleLength - 1);
    }
    setInit(extension);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_
