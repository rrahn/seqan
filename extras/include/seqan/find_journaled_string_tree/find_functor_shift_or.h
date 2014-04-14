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

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_SHIFT_OR_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_SHIFT_OR_H_

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
struct ShiftOrFunctor_;

template <typename TContainer, typename TNeedle, typename TSpec>
struct ShiftOrFunctor_<Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> >
{
    typedef Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> TFinder;
    typedef typename CallerState<TFinder>::Type TState;
    typedef unsigned int TWord;


    TWord mask;
    TState _state;
    bool _isSmallNeedle;

    ShiftOrFunctor_()
    {}

    ShiftOrFunctor_(Pattern<TNeedle, ShiftOr> & pattern)
    {
        _init(*this, pattern);
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
// Metafunction PatternSpecificTraversalSpec                         [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct PatternSpecificTraversalSpec<Pattern<TNeedle, ShiftOr> >
{
    typedef TraverserSpec<ContextPositionRight, ContextBeginLeft> Type;
};

// ----------------------------------------------------------------------------
// Metafunction CallerState
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec>
struct CallerState<Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> >
{
    typedef Pattern<TNeedle, ShiftOr> Type;
};

template <typename TContainer, typename TNeedle, typename TSpec>
struct CallerState<Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> const>
{
    typedef Pattern<TNeedle, ShiftOr> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FinderFunctor                                         [ShiftOr]
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec>
struct FinderFunctor<Finder2<THaystack, Pattern<TNeedle, ShiftOr>, TSpec> >
{
    typedef Finder2<THaystack, Pattern<TNeedle, ShiftOr>, TSpec> TFinder;
    typedef ShiftOrFunctor_<TFinder> Type;
};

template <typename THaystack, typename TNeedle, typename TSpec>
struct FinderFunctor<Finder2<THaystack, Pattern<TNeedle, ShiftOr>, TSpec> const>
{
    typedef Finder2<THaystack const, Pattern<TNeedle, ShiftOr>, TSpec> TFinder;
    typedef ShiftOrFunctor_<TFinder> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getCallerState()                                         [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec>
inline typename CallerState<Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> >::Type &
getCallerState(Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> & finder)
{
    return finder._finderFunctor._state;
}

template <typename TContainer, typename TNeedle, typename TSpec>
inline typename CallerState<Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> const>::Type &
getCallerState(Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> const & finder)
{
    return finder._finderFunctor._pattern;
}

// ----------------------------------------------------------------------------
// Function setCallerState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec>
inline void
setCallerState(Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> & finder,
               typename CallerState<Finder2<TContainer, Pattern<TNeedle, ShiftOr>, TSpec> >::Type const & state)
{
    finder._finderFunctor._state = state;
}

template <typename TContainer, typename TNeedle, typename TSpec, typename TIter>
inline typename ComputeState<TContainer>::Type
compute(Finder2<TContainer, Pattern<TNeedle, ShiftOr>, DataParallel<TSpec> > & finder,
        TIter const & hystkIt)
{
    typedef Finder2<TContainer, Pattern<TNeedle, ShiftOr>, DataParallel<TSpec> > TFinder;
    typedef typename GetTraverserForFinder_<TFinder>::Type TTraverser;
    typedef typename ComputeState<TTraverser>::Type TComputeState;

    TComputeState state(false, 1);
    if (finder._finderFunctor._isSmallNeedle)
        finder._finderFunctor(state, hystkIt, BitAlgorithmSmallNeedle());
    else
        finder._finderFunctor(state, hystkIt, BitAlgorithmLongNeedle());

    return state;
}

// ----------------------------------------------------------------------------
// Function _init()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TNeedle>
inline void
_init(ShiftOrFunctor_<TFinder> & shiftOrFunctor,
      Pattern<TNeedle, ShiftOr> & pattern)
{
    typedef ShiftOrFunctor_<TFinder> TShiftOrFunctor;
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

#endif // EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_SHIFT_OR_H_
