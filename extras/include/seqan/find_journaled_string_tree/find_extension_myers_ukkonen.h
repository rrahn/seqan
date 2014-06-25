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

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_MYERS_UKKONEN_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_MYERS_UKKONEN_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// TODO(rmaerker): @weese - can this be moved to find_myers_ukkonen.h?
template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct Spec<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
    typedef TSpec Type;
};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct Spec<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const> :
    Spec<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >{};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TPattern_>
class FinderExtensionPoint<TPattern_, MyersBitVector> : public FinderExtensionPointBase<TPattern_>
{
public:
    typedef FinderExtensionPointBase<TPattern_> TSuper;

    typedef typename GetState<FinderExtensionPoint>::Type TPatternState;
    typedef typename Needle<TPattern_>::Type TNeedle;
    typedef typename TPattern_::TLargePattern TLargePattern;
    typedef typename TPatternState::TLargeState TLargeState;
    typedef typename TPattern_::TWord TWord;
    typedef typename Spec<TPattern>::Type TSpec;

    TPatternState   _state;

    TWord X, D0, HN, HP, temp;
    TWord lastBit;
    TWord carryD0, carryHP, carryHN;  // Needed for large pattern only.
    unsigned shift, limit, currentBlock;  // Needed for large pattern only.

    bool _isSmallPattern;

    template <typename TMinScore>
    FinderExtensionPoint(TPattern & pattern, TMinScore const & minScore) : TSuper(pattern)
    {
        init(*this, minScore);
    }

    template <typename TResult, typename TIterator>
    inline void
    operator()(TResult & result, TIterator const & iter, BitAlgorithmSmallNeedle const & /*tag*/)
    {
        // computing the block
        X = this->_pattern.bitMasks[ordValue((typename Value<TNeedle>::Type) getValue(iter))] | _state.VN0;

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

        TLargePattern &largePattern = *(this->_pattern.largePattern);
        TLargeState &largeState = *_state.largeState;

//        while (position(finder) < haystack_length)
//        {
        carryD0 = carryHN = 0;
        carryHP = (int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Noting with TSpec

        // if the active cell is the last of it's block, one additional block has to be calculated
        limit = largeState.lastBlock + (unsigned)(largeState.scoreMask >> (*this->_pattern.MACHINE_WORD_SIZE - 1));

        if (limit == largePattern.blockCount)
            limit--;

        shift = largePattern.blockCount * ordValue((typename Value< TNeedle >::Type) *iter);

        // computing the necessary blocks, carries between blocks following one another are stored
        for (currentBlock = 0; currentBlock <= limit; currentBlock++)
        {
            X = *this->_pattern.bitMasks[shift + currentBlock] | largeState.VN[currentBlock];

            temp = largeState.VP[currentBlock] + (X & largeState.VP[currentBlock]) + carryD0;
            if (carryD0 != (TWord)0)
                carryD0 = temp <= largeState.VP[currentBlock];
            else
                carryD0 = temp < largeState.VP[currentBlock];

            D0 = (temp ^ largeState.VP[currentBlock]) | X;
            HN = largeState.VP[currentBlock] & D0;
            HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | D0);

            X = (HP << 1) | carryHP;
            carryHP = HP >> (this->_pattern.MACHINE_WORD_SIZE - 1);

            largeState.VN[currentBlock] = X & D0;

            temp = (HN << 1) | carryHN;
            carryHN = HN >> (this->_pattern.MACHINE_WORD_SIZE - 1);

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
                largeState.scoreMask = (TWord)1 << (this->_pattern.MACHINE_WORD_SIZE - 1);
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

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition                                 [Myers]
// ----------------------------------------------------------------------------

template <typename TPattern_>
struct ContextIteratorPosition<FinderExtensionPoint<TPattern_, MyersBitVector> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                      [Myers]
// ----------------------------------------------------------------------------

template <typename TPattern_>
struct RequireFullContext<FinderExtensionPoint<TPattern_, MyersBitVector> > : False{};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                                [Myers]
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct RegisteredExtensionPoint<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
    typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern_;
    typedef FinderExtensionPoint<TPattern_, MyersBitVector> Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TPattern_>
struct GetState<FinderExtensionPoint<TPattern_, MyersBitVector> >
{
    typedef typename PatternState<TPattern_>::Type Type;
};

template <typename TPattern_>
struct GetState<FinderExtensionPoint<TPattern_, MyersBitVector> const>
{
    typedef typename PatternState<TPattern_>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function execute()
// ----------------------------------------------------------------------------

template <typename TResult, typename TPattern_, typename TContextIter>
inline void
execute(TResult & res,
        FinderExtensionPoint<TPattern_, MyersBitVector> & extensionFunctor,
        TContextIter & contextIter)
{
    if (extensionFunctor._isSmallPattern)
        extensionFunctor(res, contextIter, BitAlgorithmSmallNeedle());
    else
        extensionFunctor(res, contextIter, BitAlgorithmLongNeedle());
}

// ----------------------------------------------------------------------------
// Function getState                                           [Myers]
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline typename GetState<FinderExtensionPoint<TPattern_, MyersBitVector> >::Type &
getState(FinderExtensionPoint<TPattern_, MyersBitVector> & extensionFunctor)
{
    return extensionFunctor._state;
}

template <typename TPattern_>
inline typename GetState<FinderExtensionPoint<TPattern_, MyersBitVector> const>::Type &
getState(FinderExtensionPoint<TPattern_, MyersBitVector> const & extensionFunctor)
{
    return extensionFunctor._state;
}

// ----------------------------------------------------------------------------
// Function setState()
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline void
setState(FinderExtensionPoint<TPattern_, MyersBitVector>  & extensionFunctor,
         typename GetState<FinderExtensionPoint<TPattern_, MyersBitVector> >::Type const & state)
{
    extensionFunctor._state = state;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TPattern_>
inline void
initState(FinderExtensionPoint<TPattern_, MyersBitVector> & extensionFunctor)
{
    extensionFunctor._state = getPattern(extensionFunctor);
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TPattern_, typename TScore>
inline void
init(FinderExtensionPoint<TPattern_, MyersBitVector> & myersFunctor,
     TScore const & scoreLimit)
{
    typedef typename TPattern_::TWord TWord;

    if (isInit(myersFunctor))
        return;

    initState(myersFunctor);
    setScoreLimit(myersFunctor._state, scoreLimit);
    _patternInit(myersFunctor._pattern, myersFunctor._state, myersFunctor);  // Initialize the pattern.
    myersFunctor._isSmallPattern = myersFunctor._pattern.largePattern == NULL;
    myersFunctor.lastBit = (TWord)1 << (myersFunctor._pattern.needleSize - 1);
    setContextSize(myersFunctor, length(host(myersFunctor._pattern)) + myersFunctor._state.maxErrors);
    setInit(myersFunctor);
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_MYERS_UKKONEN_H_
