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
// The navigator for the full score dp-matrix. We need two iterators over the
// current column and the previous column. We also store the three neighboring
// cells needed for the recursion formula.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPMatrixNavigator                          [FullDPMatrix, ScoreMatrix]
// ----------------------------------------------------------------------------

// The navigator for the score matrix.
//
// This navigator runs on a FullDPMatrix while it navigates column wise.
template <typename TValue, typename THost, typename TNavigationSpec>
class DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, TNavigationSpec>
{
public:
    typedef  DPMatrix_<TValue, FullDPMatrix, THost> TDPMatrix_;
    typedef typename Pointer_<TDPMatrix_>::Type TDPMatrixPointer_;
    typedef typename Iterator<TDPMatrix_, Standard>::Type TDPMatrixIterator;

    template <typename TBandSpec,
              std::enable_if_t<std::is_same<TBandSpec, BandOff>::value, int> = 0>
    DPMatrixNavigator_(DPMatrix_<TValue, FullDPMatrix, THost> & matrix,
                       DPBandConfig<TBandSpec> const & /*band*/)
    {
        _ptrDataContainer = &matrix;
        _prevColIteratorOffset = _dataFactors(matrix)[DPMatrixDimension_::HORIZONTAL];
        _activeColIterator = begin(matrix, Standard());
        _prevColIterator = _activeColIterator - _prevColIteratorOffset;
        _laneLeap = 1;
        *_activeColIterator = TValue{};
    }

    template <typename TBandSpec,
              std::enable_if_t<!std::is_same<TBandSpec, BandOff>::value, int> = 0>
    DPMatrixNavigator_(DPMatrix_<TValue, FullDPMatrix, THost> & matrix,
                       DPBandConfig<TBandSpec> const & band)
    {
        using TMatrixSize = typename Size<DPMatrix_<TValue, SparseDPMatrix>>::Type;
        using TSignedSize = std::make_signed_t<TMatrixSize>;

        _ptrDataContainer = &matrix;
        _prevColIteratorOffset = _dataFactors(matrix)[DPMatrixDimension_::HORIZONTAL];
        // Band begins within the first row.
        if (lowerDiagonal(band) >= 0)
        {
            _laneLeap = _min(length(matrix, DPMatrixDimension_::VERTICAL), bandSize(band));
            _activeColIterator = begin(matrix, Standard()) + _dataLengths(matrix)[DPMatrixDimension_::VERTICAL] - 1;
        }
        else if (upperDiagonal(band) <= 0)  // Band begins within the first column.
        {
            _laneLeap = 1;
            _activeColIterator = begin(matrix, Standard());
        }
        else  // Band intersects with the point of origin.
        {
            TMatrixSize lengthVertical = length(matrix, DPMatrixDimension_::VERTICAL);
            int lastPos = _max(-static_cast<TSignedSize>(lengthVertical - 1), lowerDiagonal(band));
            _laneLeap = lengthVertical + lastPos;
            _activeColIterator = begin(matrix, Standard()) + _laneLeap - 1;
        }
        // Set previous iterator to same position, one column left.
        _prevColIterator = _activeColIterator - _prevColIteratorOffset;
        *_activeColIterator =TValue{};
    }

    TDPMatrixPointer_ _ptrDataContainer     = nullptr;  // Pointer to the matrix this navigator is working on.
    TValue* _prevCellDiagonal               = nullptr;  // The previous diagonal cell
    TValue* _prevCellHorizontal             = nullptr;  // The previous Horizontal cell
    TValue* _prevCellVertical               = nullptr;  // The previous Vertical cell
    int _laneLeap                           = 0;  // Stores the jump to the next column
    size_t _prevColIteratorOffset           = 0;  // Offset to reset the previous column iterator when going to the next cell.
    TDPMatrixIterator _activeColIterator    = TDPMatrixIterator();  // The active column iterator.
    TDPMatrixIterator _prevColIterator      = TDPMatrixIterator();  // The previous column iterator.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _goNextCell                                     [banded, FirstCell]
// ----------------------------------------------------------------------------

// Needed to avoid ambigious overload.
template <typename TValue, typename TDPMatrixSpec, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnTop> const &,
            FirstCell const &)
{
    // no-op
}

// Needed to avoid ambigious overload.
template <typename TValue, typename TDPMatrixSpec, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            FirstCell const &)
{
    // no-op
}

// specialized for initialization column and all other column locations.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            FirstCell const &)
{
    // no-op
}

//  specialized for all other column types located at the top.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnTop> const &,
            FirstCell const &)
{
    --dpNavigator._laneLeap;
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator = dpNavigator._activeColIterator - dpNavigator._prevColIteratorOffset;
    _scoreOfCell(*dpNavigator._prevCellHorizontal) = _scoreOfCell(*(++dpNavigator._prevColIterator));
}

template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            FirstCell const &)
{
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator = dpNavigator._activeColIterator - dpNavigator._prevColIteratorOffset;
    _scoreOfCell(*dpNavigator._prevCellHorizontal) = _scoreOfCell(*dpNavigator._prevColIterator);
}

// version for all other column types and locations.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            FirstCell const &)
{
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator = dpNavigator._activeColIterator - dpNavigator._prevColIteratorOffset;
    _scoreOfCell(*dpNavigator._prevCellDiagonal) = _scoreOfCell(*dpNavigator._prevColIterator);
    _scoreOfCell(*dpNavigator._prevCellHorizontal) = _scoreOfCell(*(++dpNavigator._prevColIterator));
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                   [unbanded, FirstCell]
// ----------------------------------------------------------------------------

// specialized for initalization column.
template <typename TValue, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            FirstCell const &)
{
    // no-op
}

// all other column types.
template <typename TValue, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            FirstCell const &)
{
    // Set to begin of column.
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator = dpNavigator._activeColIterator - dpNavigator._prevColIteratorOffset;
    // Cache prevDiagH value. Becomes the next diagonal reference value.
    _scoreOfCell(*dpNavigator._prevCellHorizontal) = _scoreOfCell(*dpNavigator._prevColIterator);

}

// ----------------------------------------------------------------------------
// Function _goNextCell                                     [banded, InnerCell]
// ----------------------------------------------------------------------------

// specilized for the initialization column and all column locations.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            InnerCell const &)
{
    _cachePrevVertical(dpNavigator);
    ++dpNavigator._activeColIterator;
}

// specialized for the standard band processing.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            InnerCell const &)
{
    _cachePrevVertical(dpNavigator);
    _scoreOfCell(*dpNavigator._prevCellDiagonal) = _scoreOfCell(*dpNavigator._prevCellHorizontal);
    ++dpNavigator._activeColIterator;
    _scoreOfCell(*dpNavigator._prevCellHorizontal) = _scoreOfCell(*(++dpNavigator._prevColIterator));
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                   [unbanded, InnerCell]
// ----------------------------------------------------------------------------

// specialized for the initialization column.
template <typename TValue, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            InnerCell const &)
{
    _cachePrevVertical(dpNavigator);
    ++dpNavigator._activeColIterator;
}

// version for the all other column types.
template <typename TValue, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            InnerCell const &)
{
    _cachePrevVertical(dpNavigator);
    // Cache prevDiag from prevDiagH;
    _scoreOfCell(*dpNavigator._prevCellDiagonal) = _scoreOfCell(*dpNavigator._prevCellHorizontal);
    // Cache prevDiagH from current;
    _scoreOfCell(*dpNavigator._prevCellHorizontal) = _scoreOfCell(*(++dpNavigator._prevColIterator));
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                      [banded, LastCell]
// ----------------------------------------------------------------------------

// specilaized for initialization column and bottom column.
template <typename TValue, typename TDPMatrixSpec, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom> const &,
            LastCell const &)
{
    // dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    _cachePrevVertical(dpNavigator);
    ++dpNavigator._activeColIterator;
}

// specilaized for initialization column and all other column locations.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            LastCell const &)
{
    _cachePrevVertical(dpNavigator);
    ++dpNavigator._activeColIterator;
}

// specialized for all column types and bottom column.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnBottom> const &,
            LastCell const &)
{
    _cachePrevVertical(dpNavigator);
    _scoreOfCell(*dpNavigator._prevCellDiagonal) = _scoreOfCell(*dpNavigator._prevCellHorizontal);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._prevColIterator;
    ++dpNavigator._laneLeap;
}

// generic case for all other column types and column locations.
template <typename TValue, typename TDPMatrixSpec, typename THost,
          typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, TDPMatrixSpec, THost>, DPScoreMatrix, NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            LastCell const &)
{
    _cachePrevVertical(dpNavigator);
    _scoreOfCell(*dpNavigator._prevCellDiagonal) = _scoreOfCell(*dpNavigator._prevCellHorizontal);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._prevColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                    [unbanded, LastCell]
// ----------------------------------------------------------------------------

// specilaized for initialization column.
template <typename TValue, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            LastCell const &)
{
    _cachePrevVertical(dpNavigator);
    ++dpNavigator._activeColIterator;
}

// version for all other column types.
template <typename TValue, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            LastCell const &)
{
    // Cache the vertical cells: prevDiagV and vertical.
    _cachePrevVertical(dpNavigator);
    // Cache prevDiag from prevDiagH
    _scoreOfCell(*dpNavigator._prevCellDiagonal) = _scoreOfCell(*dpNavigator._prevCellHorizontal);
    ++dpNavigator._activeColIterator; // go to next cell.
    ++dpNavigator._prevColIterator; // go to next cell.
}

// ----------------------------------------------------------------------------
// Function _cachePrevVertical()
// ----------------------------------------------------------------------------

template <typename TValue, typename TGapSpec, typename TMatrixSpec, typename THost, typename TNavigationSpec>
inline void
_cachePrevVertical(DPMatrixNavigator_<DPMatrix_<DPCell_<TValue, TGapSpec>, TMatrixSpec, THost>, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    // Cache prevVerical.
    *dpNavigator._prevCellVertical = *dpNavigator._activeColIterator;
}

template <typename TValue, typename TMatrixSpec, typename THost, typename TNavigationSpec>
inline void
_cachePrevVertical(DPMatrixNavigator_<DPMatrix_<DPCell_<TValue, AffineGaps>, TMatrixSpec, THost>, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    // Cache prevDiagV.
    _verticalScoreOfCell(*dpNavigator._prevCellVertical) = _verticalScoreOfCell(*dpNavigator._activeColIterator);
    // Cache prevVerical.
    _scoreOfCell(*dpNavigator._prevCellVertical) = _scoreOfCell(*dpNavigator._activeColIterator);
}

// ----------------------------------------------------------------------------
// Function previousCellDiagonal()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> >::Type
previousCellDiagonal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    return *dpNavigator._prevCellDiagonal;
}

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const>::Type
previousCellDiagonal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const & dpNavigator)
{
    return *dpNavigator._prevCellDiagonal;
}

// ----------------------------------------------------------------------------
// Function previousCellHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> >::Type
previousCellHorizontal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    return *dpNavigator._prevColIterator;
}

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const>::Type
previousCellHorizontal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const & dpNavigator)
{
    return *dpNavigator._prevColIterator;
}

// ----------------------------------------------------------------------------
// Function previousCellVertical()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> >::Type
previousCellVertical(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    return *dpNavigator._prevCellVertical;
}

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const>::Type
previousCellVertical(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const & dpNavigator)
{
    return *dpNavigator._prevCellVertical;
}

// ----------------------------------------------------------------------------
// Function setCachedCells()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename TNavigationSpec>
inline void
setCachedCells(DPMatrixNavigator_<DPMatrix_<TValue, TMatrixSpec>, DPScoreMatrix, TNavigationSpec> & dpNavigator,
               TValue & prevDiag,
               TValue & prevHor,
               TValue & prevVer)
{
    dpNavigator._prevCellDiagonal = &prevDiag;
    dpNavigator._prevCellHorizontal = &prevHor;
    dpNavigator._prevCellVertical = &prevVer;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_H_
