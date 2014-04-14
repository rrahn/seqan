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
// Implements the meyers ukkonen algorithm for approximate pattern matching
// with edit distance.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_MYERS_UKKONEN_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_MYERS_UKKONEN_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename T>
struct MyersUkkonenFunctor_;

template <typename TContainer, typename TNeedle, typename TSpec, typename THostSpec, typename TPatternBeginSpec, typename TFinderSpec>
struct MyersUkkonenFunctor_<Finder2<TContainer, Pattern<TNeedle, Myers<TSpec, THostSpec, TPatternBeginSpec> >, DataParallel<TFinderSpec> > >
{
    typedef Pattern<TNeedle, Myers<TSpec, THostSpec, TPatternBeginSpec> > TPattern;
    typedef typename PatternState<TPattern>::Type TPatternState;
    typedef Finder2<TContainer, TPattern, DataParallel<TFinderSpec> > TFinder;

    typedef MyersLargePattern_<TNeedle, TSpec> TLargePattern;
    typedef MyersLargeState_<TNeedle, TSpec> TLargeState;
    typedef typename TPattern::TWord TWord;

    TPattern        _pattern;
    TPatternState   _state;

    TWord X, D0, HN, HP, temp;
    TWord lastBit;
    TWord carryD0, carryHP, carryHN;  // Needed for large pattern only.
    unsigned shift, limit, currentBlock;  // Needed for large pattern only.

    bool _isSmallPattern;

    MyersUkkonenFunctor_()
    {}

    template <typename TErrors>
    MyersUkkonenFunctor_(TPattern & pattern, TErrors const & errors)
    {
        _init(*this, pattern, errors);
    }

    template <typename TResult, typename TIterator>
    inline void
    operator()(TResult & result, TIterator const & iter, BitAlgorithmSmallNeedle const & /*tag*/)
    {
        // computing the block
        X = _pattern.bitMasks[ordValue((typename Value<TNeedle>::Type) getValue(iter))] | _state.VN0;

        D0 = ((_state.VP0 + (X & _state.VP0)) ^ _state.VP0) | X;
        HN = _state.VP0 & D0;
        HP = _state.VN0 | ~(_state.VP0 | D0);
        X = (HP << 1) | (TWord)(int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Nothing by TSpec
        _state.VN0 = X & D0;
        _state.VP0 = (HN << 1) | ~(X | D0);

        if ((HP & lastBit) != (TWord)0)
            _state.errors++;
        else if ((HN & lastBit) != (TWord)0)
            _state.errors--;

        result.i1 = _state.errors <= _state.maxErrors;
    }

    template <typename TResult, typename TIterator>
    inline void
    operator()(TResult & result, TIterator const & iter, BitAlgorithmLongNeedle const & /*tag*/)
    {
        // computing the block

        TLargePattern &largePattern = *_pattern.largePattern;
        TLargeState &largeState = *_state.largeState;

//        while (position(finder) < haystack_length)
//        {
        carryD0 = carryHN = 0;
        carryHP = (int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Noting with TSpec

        // if the active cell is the last of it's block, one additional block has to be calculated
        limit = largeState.lastBlock + (unsigned)(largeState.scoreMask >> (_pattern.MACHINE_WORD_SIZE - 1));

        if (limit == largePattern.blockCount)
            limit--;

        shift = largePattern.blockCount * ordValue((typename Value< TNeedle >::Type) *iter);

        // computing the necessary blocks, carries between blocks following one another are stored
        for (currentBlock = 0; currentBlock <= limit; currentBlock++)
        {
            X = _pattern.bitMasks[shift + currentBlock] | largeState.VN[currentBlock];

            temp = largeState.VP[currentBlock] + (X & largeState.VP[currentBlock]) + carryD0;
            if (carryD0 != (TWord)0)
                carryD0 = temp <= largeState.VP[currentBlock];
            else
                carryD0 = temp < largeState.VP[currentBlock];

            D0 = (temp ^ largeState.VP[currentBlock]) | X;
            HN = largeState.VP[currentBlock] & D0;
            HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | D0);

            X = (HP << 1) | carryHP;
            carryHP = HP >> (_pattern.MACHINE_WORD_SIZE - 1);

            largeState.VN[currentBlock] = X & D0;

            temp = (HN << 1) | carryHN;
            carryHN = HN >> (_pattern.MACHINE_WORD_SIZE - 1);

            largeState.VP[currentBlock] = temp | ~(X | D0);

            // if the current block is the one containing the last active cell
            // the new number of errors is computed
            if (currentBlock == largeState.lastBlock) {
                if ((HP & largeState.scoreMask) != (TWord)0)
                    _state.errors++;
                else if ((HN & largeState.scoreMask) != (TWord)0)
                    _state.errors--;
            }
        }

        // updating the last active cell
        while (_state.errors > _state.maxErrors) {
            if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                _state.errors--;
            else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                _state.errors++;

            largeState.scoreMask >>= 1;
            if (largeState.scoreMask == (TWord)0)
            {
                largeState.lastBlock--;
                if (IsSameType<TSpec, FindPrefix>::VALUE && largeState.lastBlock == (unsigned)-1)
                    break;
                largeState.scoreMask = (TWord)1 << (_pattern.MACHINE_WORD_SIZE - 1);
            }
        }

        if ((largeState.scoreMask == largePattern.finalScoreMask) && (largeState.lastBlock == largePattern.blockCount - 1))
        {
//                _setFinderEnd(finder);
//                if (IsSameType<TSpec, FindPrefix>::VALUE)
//                {
//                    _setFinderLength(finder, endPosition(finder));
//                }
            result.i1 = true;
//                return true;
        }
        else {
            largeState.scoreMask <<= 1;
            if (!largeState.scoreMask) {
                largeState.scoreMask = 1;
                largeState.lastBlock++;
            }

            if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                _state.errors++;
            else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                _state.errors--;
        }

    //      SEQAN_ASSERT (state.errors >= 0);
//
//            gonext(finder);
    }

//        return false;
//    }

};


// ============================================================================
// Metafunctions
// ============================================================================

template <typename TFinder>
struct ErrorsSupported<MyersUkkonenFunctor_<TFinder> > : True{};

// ----------------------------------------------------------------------------
// Metafunction PatternSpecificTraversalSpec                         [ShiftAnd]
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename THasState, typename TBeginPatternSpec>
struct PatternSpecificTraversalSpec<Pattern<TNeedle, Myers<TSpec, THasState, TBeginPatternSpec> > >
{
    typedef TraverserSpec<ContextPositionRight, ContextBeginLeft> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FinderFunctor
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
struct FinderFunctor<Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> >
{
    typedef Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> TFinder_;
    typedef MyersUkkonenFunctor_<TFinder_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction InitStateForFinder
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
struct InitStateForFinder<Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, DataParallel<TFinderSpec> > >
{
    typedef FinderInitializationState<Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> > > Type;
};

// ----------------------------------------------------------------------------
// Metafunction CallerState
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
struct CallerState<Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> >
{
    typedef Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> > TPattern_;
    typedef typename PatternState<TPattern_>::Type Type;
};

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
struct CallerState<Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> const>
{
    typedef Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> > TPattern_;
    typedef typename PatternState<TPattern_>::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function compute()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TPatternSpec, typename THasState,
          typename TFindBeginPatternSpec, typename TSpec, typename TIterator>
inline typename ComputeState<TContainer>::Type
compute(Finder2<TContainer, Pattern<TNeedle, Myers<TPatternSpec, THasState, TFindBeginPatternSpec> >, DataParallel<TSpec> > & finder,
        TIterator const & iter)
{
    typedef Pattern<TNeedle, Myers<TPatternSpec, THasState, TFindBeginPatternSpec> > TPattern;
    typedef typename PatternSpecificTraversalSpec<TPattern>::Type TTraversalSpec;
    typedef Traverser<TContainer, TTraversalSpec> TTraverser;

    typedef typename ComputeState<TTraverser>::Type TComputeState;

    TComputeState state(false, 1);
    if (finder._finderFunctor._isSmallPattern)
        finder._finderFunctor(state, iter, BitAlgorithmSmallNeedle());
    else
        finder._finderFunctor(state, iter, BitAlgorithmLongNeedle());
    return state;
}

// ----------------------------------------------------------------------------
// Function getCallerState                                              [Myers]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
inline typename CallerState<Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> >::Type &
getCallerState(Finder2<TContainer, Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec > & finder)
{
    return finder._finderFunctor._state;
}

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
inline typename CallerState<Finder2<TContainer,  Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> const>::Type &
getCallerState(Finder2<TContainer, Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> const & finder)
{
    return finder._finderFunctor._state;
}

// ----------------------------------------------------------------------------
// Function setCallerState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec, typename TFinderSpec>
inline void
setCallerState(Finder2<TContainer, Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> & finder,
               typename CallerState<Finder2<TContainer, Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, TFinderSpec> >::Type const & state)
{
    finder._finderFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function setScoreLimit()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec,
          typename TFinderSpec, typename TScore>
inline void
setScoreLimit(MyersUkkonenFunctor_<Finder2<TContainer, Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> >, DataParallel<TFinderSpec> > > & myersFunctor,
              TScore const & score)
{
    setScoreLimit(myersFunctor._state, score);
}

// ----------------------------------------------------------------------------
// Function _reinit()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TNeedle, typename TSpec, typename THasHost, typename TPatternBeginSpec,
          typename TScore>
inline void
_init(MyersUkkonenFunctor_<TFinder> & myersFunctor,
      Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> > & pattern,
      TScore const & scoreLimit)
{
    typedef Pattern<TNeedle, Myers<TSpec, THasHost, TPatternBeginSpec> > TPattern;
    typedef typename TPattern::TWord TWord;

    myersFunctor._pattern = pattern;
    myersFunctor._state = pattern;

    setScoreLimit(myersFunctor, scoreLimit);
    _patternInit(myersFunctor._pattern, myersFunctor._state, myersFunctor);  // Initialize the pattern.
    myersFunctor._isSmallPattern = myersFunctor._pattern.largePattern == NULL;
    myersFunctor.lastBit = (TWord)1 << (pattern.needleSize - 1);
}


}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_MYERS_UKKONEN_H_
