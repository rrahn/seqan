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
// Defines the methods to compute the score when using linear gap costs.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_

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
// Function _doComputeScore                 [RecursionDirectionAll, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, 
          typename TSequenceHValue, 
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TRecursionDirection,
          typename TAlgorithm, 
          typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                TRecursionDirection const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;

    TScoreValue sv = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);

    TTraceValue tv = _dpMax((_scoreOfCell(activeCell) = -13436436), // TODO:: replace by infinity
                            sv,
                            +TraceBitMap_<TScoreValue>::NONE,
                            TraceBitMap_<TScoreValue>::HORIZONTAL | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                            TTracebackConfig(),
                            typename _activate<std::decay_t<TRecursionDirection>, DPHorizontal>::Type());
    
    tv = _dpMax(_scoreOfCell(activeCell),
                static_cast<TScoreValue>(_scoreOfCell(previousVertical) +
                                         scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal)),
                tv,
                TraceBitMap_<TScoreValue>::VERTICAL | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                TTracebackConfig(),
                typename _activate<std::decay_t<TRecursionDirection>, DPVertical>::Type());

    return _dpMax(_scoreOfCell(activeCell),
                  static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                           score(scoringScheme, seqHVal, seqVVal)),
                  tv,
                  +TraceBitMap_<TScoreValue>::DIAGONAL,
                  TTracebackConfig(),
                  typename _activate<std::decay_t<TRecursionDirection>, DPDiagonal>::Type());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
