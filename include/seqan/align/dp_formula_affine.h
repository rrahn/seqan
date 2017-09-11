// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Implements the affine gap cost functions.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _internalComputeScore          [Vertical vs Horizontal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue>>>,
                            typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOff const &)
{
    using std::max;
    _scoreOfCell(activeCell) = max(_horizontalScoreOfCell(activeCell), _verticalScoreOfCell(activeCell));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue>>,
                            typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOff const &)
{
    _scoreOfCell(activeCell) = max(_horizontalScoreOfCell(activeCell), _verticalScoreOfCell(activeCell));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    return (activeCell._verticalScore < activeCell._horizontalScore)
        ? (activeCell._score = activeCell._horizontalScore, TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX)
        : (activeCell._score = activeCell._verticalScore, TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    TScoreValue cmp = cmpGt(activeCell._horizontalScore, activeCell._verticalScore);
    activeCell._score = blend(activeCell._verticalScore, activeCell._horizontalScore, cmp);
    return blend(TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                 TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                 cmp);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    return (activeCell._horizontalScore <= activeCell._verticalScore)
        ? ((activeCell._horizontalScore == activeCell._verticalScore)
            ? (activeCell._score = activeCell._horizontalScore,
               TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX |
               TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX)
            : (activeCell._score = activeCell._verticalScore, TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX))
        : (activeCell._score = activeCell._horizontalScore, TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    TScoreValue cmpG = cmpGt(activeCell._horizontalScore, activeCell._verticalScore);
    TScoreValue cmpE = cmpEq(activeCell._horizontalScore, activeCell._verticalScore);
    activeCell._score = blend(activeCell._verticalScore, activeCell._horizontalScore, cmpG);

    return blend(blend(TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                       TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                       cmpG),
                 TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                 cmpE);
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                 [RecursionAllDirection, AffineGaps]
// ----------------------------------------------------------------------------

// if veritcal or horizontal recursion direction is considered,
// the score must be assigned and the trace can be used
template <typename TScoreValue,
          typename TTraceValue>
inline typename TraceBitMap_<TScoreValue>::Type
_assignScoreHV(TScoreValue & score, TScoreValue const & newScore, TTraceValue const & trace, True)
{
    score = newScore;
    return trace;
}

// if not, the score must be set to infinity and the trace must not be [DIRECTION]_OPEN
// but none.
template <typename TScoreValue,
          typename TTraceValue>
inline typename TraceBitMap_<TScoreValue>::Type
_assignScoreHV(TScoreValue & score, TScoreValue const &, TTraceValue const &, False)
{
    score = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    return +TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue,
          typename TSequenceHValue, 
          typename TSequenceVValue, 
          typename TScoringScheme,
          typename TRecursionDirection,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                TRecursionDirection const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;

    // Compute horizontal direction.
    // this needs to be computed first because _horizontalScoreOfCell returns *activeColIterator in a SparseMatrix
    TScoreValue sv = _horizontalScoreOfCell(previousHorizontal) +
                     scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _dpMax(_horizontalScoreOfCell(activeCell),
                            sv,
                            _assignScoreHV(_horizontalScoreOfCell(activeCell),
                                           static_cast<TScoreValue>(_scoreOfCell(previousHorizontal) +
                                                                    scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal)),
                                           +TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                                           typename _activate<std::decay_t<TRecursionDirection>, DPHorizontal>::Type()),
                            +TraceBitMap_<TScoreValue>::HORIZONTAL,
                            TTracebackConfig(),
                            typename _activate<std::decay_t<TRecursionDirection>, DPHorizontal>::Type());

    // Compute vertical direction.
    tv |= _dpMax(_verticalScoreOfCell(activeCell),
                 static_cast<TScoreValue>(_verticalScoreOfCell(previousVertical) +
                                          scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal)),
                 _assignScoreHV(_verticalScoreOfCell(activeCell),
                                static_cast<TScoreValue>(_scoreOfCell(previousVertical) +
                                                         scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal)),
                                +TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                                typename _activate<std::decay_t<TRecursionDirection>, DPVertical>::Type()),
                 +TraceBitMap_<TScoreValue>::VERTICAL,
                 TTracebackConfig(),
                 typename _activate<std::decay_t<TRecursionDirection>, DPVertical>::Type());

    // Get max from horiztonal and/or vertical direction
    TTraceValue tvMax = _internalComputeScore(activeCell, TTracebackConfig());

    // and compare with diagonal direction.
    return _dpMax(_scoreOfCell(activeCell),
                  static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                           score(scoringScheme, seqHVal, seqVVal)),
                  tvMax | tv,
                  TraceBitMap_<TScoreValue>::DIAGONAL | tv, 
                  TTracebackConfig(),
                  typename _activate<std::decay_t<TRecursionDirection>, DPDiagonal>::Type());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
