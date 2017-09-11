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
// Defines the recursion formula for the dp-alignment algorithms.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag RecursionDirectionDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionDiagonal_;
typedef Tag<RecursionDirectionDiagonal_> RecursionDirectionDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionHorizontal
// ----------------------------------------------------------------------------

struct RecursionDirectionHorizontal_;
typedef Tag<RecursionDirectionHorizontal_> RecursionDirectionHorizontal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionVertical
// ----------------------------------------------------------------------------

struct RecursionDirectionVertical_;
typedef Tag<RecursionDirectionVertical_> RecursionDirectionVertical;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionAll
// ----------------------------------------------------------------------------

struct RecursionDirectionAll_;
typedef Tag<RecursionDirectionAll_> RecursionDirectionAll;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionUpperDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionUpperDiagonal_;
typedef Tag<RecursionDirectionUpperDiagonal_> RecursionDirectionUpperDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionLowerDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionLowerDiagonal_;
typedef Tag<RecursionDirectionLowerDiagonal_> RecursionDirectionLowerDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionZero
// ----------------------------------------------------------------------------

struct RecursionDirectionZero_;
typedef Tag<RecursionDirectionZero_> RecursionDirectionZero;

// ============================================================================
// Metafunctions
// ============================================================================

// Helper typedef to get the correct score value type from the score-matrix navigator.
template <typename TCellTuple>
using ExtractedScoreValueType_ = std::decay_t<decltype(_scoreOfCell(std::get<0>(std::declval<TCellTuple>())))>;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _computeScore
// ----------------------------------------------------------------------------

template <typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TRecursionDirection,
          typename TDPProfile>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<ExtractedScoreValueType_<TRecursionCellTuple>>>>,
                            typename TraceBitMap_<ExtractedScoreValueType_<TRecursionCellTuple>>::Type)
_computeScore(TRecursionCellTuple && recursionCells,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              TRecursionDirection const & recDir,
              TDPProfile const & dpProfile)
{
    auto traceDir = _doComputeScore(std::get<0>(recursionCells),
                                    std::get<1>(recursionCells),
                                    std::get<2>(recursionCells),
                                    std::get<3>(recursionCells),
                                    seqHVal, seqVVal, scoringScheme, recDir, dpProfile);
    if (IsLocalAlignment_<TDPProfile>::VALUE)
    {
        if (_scoreOfCell(std::get<0>(recursionCells)) <= 0)
        {
            _setScoreOfCell(std::get<0>(recursionCells), static_cast<ExtractedScoreValueType_<TRecursionCellTuple>>(0));
            return TraceBitMap_<ExtractedScoreValueType_<TRecursionCellTuple>>::NONE;
        }
    }
    return traceDir;
}

template <typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TRecursionDirection,
          typename TDPProfile>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<ExtractedScoreValueType_<TRecursionCellTuple>>>,
                            typename TraceBitMap_<ExtractedScoreValueType_<TRecursionCellTuple>>::Type)
_computeScore(TRecursionCellTuple && recursionCells,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              TRecursionDirection const & recDir,
              TDPProfile const & dpProfile)
{
    using TScoreValue = ExtractedScoreValueType_<TRecursionCellTuple>;

    auto traceDir = _doComputeScore(std::get<0>(recursionCells),
                                    std::get<1>(recursionCells),
                                    std::get<2>(recursionCells),
                                    std::get<3>(recursionCells),
                                    seqHVal, seqVVal, scoringScheme, recDir, dpProfile);
    if (IsLocalAlignment_<TDPProfile>::VALUE)
    {
        TScoreValue cmp = cmpGt(createVector<TScoreValue>(1), _scoreOfCell(std::get<0>(recursionCells)));
        _setScoreOfCell(std::get<0>(recursionCells), TraceBitMap_<TScoreValue>::NONE, cmp);
        return blend(traceDir, TraceBitMap_<TScoreValue>::NONE, cmp);
    }
    return traceDir;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionDirectionDiagonal]
// ----------------------------------------------------------------------------
// Independent of gap cost model.
template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline auto
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, AffineGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionDiagonal const &,
                DPProfile_<TAlgoTag, AffineGaps, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), False(), createVector<TScoreValue>(-1));

    if (!IsTracebackEnabled_<DPProfile_<TAlgoTag, AffineGaps, TTraceFlag>>::VALUE)
        return TraceBitMap_<TScoreValue>::NONE;
    return TraceBitMap_<TScoreValue>::DIAGONAL;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                            [RecursionDirectionZero]
// ----------------------------------------------------------------------------
// Independent of gap cost model.
template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline auto
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, AffineGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, AffineGaps> const & /*previousVertical*/,
                TSequenceHValue const & /*seqHVal*/,
                TSequenceVValue const & /*seqVVal*/,
                TScoringScheme const & /*scoringScheme*/,
                RecursionDirectionZero const &,
                DPProfile_<TAlgoTag, AffineGaps, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = createVector<TScoreValue>(0);
    return TraceBitMap_<TScoreValue>::NONE;
}
// Independent of gap cost model.
template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TDPProfile>
inline auto
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & previousDiagonal,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionDiagonal const &,
                TDPProfile const &)
{
    _scoreOfCell(activeCell) = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), False(), createVector<TScoreValue>(-1));

    if (!IsTracebackEnabled_<TDPProfile>::VALUE)
        return TraceBitMap_<TScoreValue>::NONE;
    return TraceBitMap_<TScoreValue>::DIAGONAL;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                            [RecursionDirectionZero]
// ----------------------------------------------------------------------------
// Independent of gap cost model.
template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline auto
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousVertical*/,
                TSequenceHValue const & /*seqHVal*/,
                TSequenceVValue const & /*seqVVal*/,
                TScoringScheme const & /*scoringScheme*/,
                RecursionDirectionZero const &,
                DPProfile_<TAlgoTag, DynamicGaps, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = createVector<TScoreValue>(0);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline auto
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, LinearGaps> const & /*previousVertical*/,
                TSequenceHValue const & /*seqHVal*/,
                TSequenceVValue const & /*seqVVal*/,
                TScoringScheme const & /*scoringScheme*/,
                RecursionDirectionZero const &,
                DPProfile_<TAlgoTag, LinearGaps, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = createVector<TScoreValue>(0);
    return TraceBitMap_<TScoreValue>::NONE;
}
// Independent of gap cost model.
template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline auto
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, LinearGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionDiagonal const &,
                DPProfile_<TAlgoTag, LinearGaps, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), False(), createVector<TScoreValue>(-1));

    if (!IsTracebackEnabled_<DPProfile_<TAlgoTag, LinearGaps, TTraceFlag>>::VALUE)
        return TraceBitMap_<TScoreValue>::NONE;
    return TraceBitMap_<TScoreValue>::DIAGONAL;
}

// Metafunction 
struct DPVertical_;
typedef Tag<DPVertical_> DPVertical;

struct DPHorizontal_;
typedef Tag<DPHorizontal_> DPHorizontal;

struct DPDiagonal_;
typedef Tag<DPDiagonal_> DPDiagonal;

struct DPHorizontalOrVertical_;
typedef Tag<DPHorizontalOrVertical_> DPHorizontalOrVertical;

template <typename TRecursionDirection, typename TMaxComputation>
struct _activate
{
    typedef False Type; 
};

template <typename TMaxComputation>
struct _activate<RecursionDirectionAll, TMaxComputation>
{
    typedef True Type; 
};

template <typename TMaxComputation>
struct _activate<RecursionDirectionDiagonal, TMaxComputation>
{
    typedef False Type; 
};

template <>
struct _activate<RecursionDirectionDiagonal, DPDiagonal>
{
    typedef True Type; 
};

template <typename TMaxComputation>
struct _activate<RecursionDirectionHorizontal, TMaxComputation>
{
    typedef False Type; 
};

template <>
struct _activate<RecursionDirectionHorizontal, DPHorizontal>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionHorizontal, DPHorizontalOrVertical>
{
    typedef True Type; 
};

template <typename TMaxComputation>
struct _activate<RecursionDirectionVertical, TMaxComputation>
{
    typedef False Type; 
};

template <>
struct _activate<RecursionDirectionVertical, DPVertical>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionVertical, DPHorizontalOrVertical>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionUpperDiagonal, DPDiagonal>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionUpperDiagonal, DPHorizontal>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionUpperDiagonal, DPVertical>
{
    typedef False Type; 
};

template <>
struct _activate<RecursionDirectionUpperDiagonal, DPHorizontalOrVertical>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionLowerDiagonal, DPDiagonal>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionLowerDiagonal, DPHorizontal>
{
    typedef False Type; 
};

template <>
struct _activate<RecursionDirectionLowerDiagonal, DPVertical>
{
    typedef True Type; 
};

template <>
struct _activate<RecursionDirectionLowerDiagonal, DPHorizontalOrVertical>
{
    typedef True Type; 
};

// ----------------------------------------------------------------------------
// Function _dpMax                        
// ----------------------------------------------------------------------------

// Brief: This function determines the max between two scores.
//
// Detail:
// Doing so, it assigns the maximal Score to the lhs value (activeScore)
// and depending on which score was maximal returns the correspong trace
// (rhsTrace if activeScore < rhsScore, lhsTrace else)

// The function is overloaded for the follwoing traceback options:
// TracebackOff  - No SIMD
// TracebackOff  - SIMD
// SingleTrace   - No SIMD
// SingleTrace   - SIMD
// CompleteTrace - No SIMD
// CompleteTrace - SIMD

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_dpMax(TScoreValue & activeScore,
       TScoreValue const & rhsScore,
       TTraceValueL,
       TTraceValueR,
       TracebackOff const &,
       True)
{
    using std::max;
    activeScore = max(activeScore, rhsScore);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_dpMax(TScoreValue & activeScore,
       TScoreValue const & rhsScore,
       TTraceValueL,
       TTraceValueR,
       TracebackOff const &,
       True)
{
    TScoreValue cmp = cmpGt(rhsScore, activeScore);
    activeScore = blend(activeScore, rhsScore, cmp);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_dpMax(TScoreValue & activeScore,
       TScoreValue const & rhsScore,
       TTraceValueL lhsTrace,
       TTraceValueR rhsTrace,
       TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
       True)
{

    return (activeScore <= rhsScore) 
    ? (activeScore = rhsScore, rhsTrace) 
    : (lhsTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_dpMax(TScoreValue & activeScore,
       TScoreValue const & rhsScore,
       TTraceValueL lhsTrace,
       TTraceValueR rhsTrace,
       TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
       True)
{
    TScoreValue cmp = cmpGt(rhsScore, activeScore);
    activeScore = blend(activeScore, rhsScore, cmp);
    return blend(lhsTrace, rhsTrace, cmp);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
    inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_dpMax(TScoreValue & activeScore,
       TScoreValue const & rhsScore,
       TTraceValueL lhsTrace,
       TTraceValueR rhsTrace,
       TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const &,
       True)
{
    return (activeScore <= rhsScore) ?
        ((activeScore == rhsScore) ? (rhsTrace | lhsTrace) :
                                     (activeScore = rhsScore, rhsTrace)) :
        (lhsTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_dpMax(TScoreValue & activeScore,
       TScoreValue const & rhsScore,
       TTraceValueL const & lhsTrace,
       TTraceValueR const & rhsTrace,
       TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const &,
       True)
{
    // Check for greater values.
    TScoreValue cmp = cmpGt(activeScore, rhsScore);  // cmp greater
    activeScore = blend(rhsScore, activeScore, cmp);  // activeScore

    // Check for equality.
    return blend(blend(rhsTrace, lhsTrace, cmp), lhsTrace | rhsTrace, cmpEq(rhsScore, activeScore));
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TTraceConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_dpMax(TScoreValue const &,
       TScoreValue const &,
       TTraceValueL const & leftTrace,
       TTraceValueR const &,
       TTraceConfig const &,
       False)
{
    return leftTrace;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
