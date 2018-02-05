// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: René Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/stream.h>

#include "test_align_parallel_algorithm_banded.h"

SEQAN_BEGIN_TESTSUITE(test_simd_vector)
{
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_construction);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_begin);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_end);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_increment_pre);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_increment_post);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_increment_offset);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_decrement_pre);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_decrement_post);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_decrement_offset);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_difference);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_comparison);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_dereferencing);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_container);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_position);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_coordinate);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_column_index);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_row_index);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_predecessor_left);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_predecessor_above);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_successor_right);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_successor_below);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_successor_right);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_successor_below);
    SEQAN_CALL_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_transform_to_grid);
    // First test some basic data structures.
// #if defined(SEQAN_SEQANSIMD_ENABLED) && defined(__SSE4_1__)
//     SEQAN_CALL_TEST(test_simd_transpose_8x8);
//     SEQAN_CALL_TEST(test_simd_transpose_16x16);
//     SEQAN_CALL_TEST(test_simd_types);
// #endif  // defined(SEQAN_SEQANSIMD_ENABLED) && defined(__SSE4_1__)

// #if defined(SEQAN_SEQANSIMD_ENABLED) && defined(__AVX2__)
//     SEQAN_CALL_TEST(test_simd_transpose_32x32);
// #endif  // defined(SEQAN_SEQANSIMD_ENABLED) && defined(__AVX2__)

    return seqan::TestSystem::runAll();
}
SEQAN_END_TESTSUITE
