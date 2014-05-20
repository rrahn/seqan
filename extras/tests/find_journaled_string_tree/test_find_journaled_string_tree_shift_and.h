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
// Implements test for data prallel finder with the shift and.
// ==========================================================================

#ifndef EXTRAS_TESTS_FIND_JOURNALED_STRING_TREE_TEST_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_
#define EXTRAS_TESTS_FIND_JOURNALED_STRING_TREE_TEST_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_

#include <seqan/basic.h>
#include <seqan/find_journaled_string_tree.h>

#include "test_find_journaled_string_tree_compare_results.h"

SEQAN_DEFINE_TEST(test_find_journaled_string_tree_ref_occ_shift_and)
{
    using namespace seqan;

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 0, -1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 1, -1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, -1, 0, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 0, -1, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 1, -1, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, -1, 10, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 0, -1, 35, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 1, -1, 35, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, -1, 35, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 0, -1, 80, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 1, -1, 80, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, -1, 80, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 0, -1, 96, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 1, -1, 96, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, -1, 96, 10));
}

SEQAN_DEFINE_TEST(test_find_journaled_string_tree_snp_occ_shift_and)
{
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 0, 0, 0, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 0, 1, 1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 0, 5, 2, 0, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 0, 0, 3, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 0, 1, 4, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 0, 5, 5, 10, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 0, 0, 32, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 1, 2, 32, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, 4, 32, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 0, 0, -1, 5, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 0, 1, -1, 33, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 0, 5, -1, 50, 10));
}

SEQAN_DEFINE_TEST(test_find_journaled_string_tree_ins_occ_shift_and)
{
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 2, 0, 0, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 2, 1, 1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 2, 5, 2, 0, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 2, 0, 3, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 2, 1, 4, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 2, 5, 5, 10, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 2, 0, 0, 32, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 2, 1, 2, 32, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 2, 5, 4, 32, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 2, 0, -1, 5, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 2, 1, -1, 33, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 2, 5, -1, 50, 10));
}

SEQAN_DEFINE_TEST(test_find_journaled_string_tree_del_occ_shift_and)
{
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 1, 0, 0, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 1, 1, 1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 1, 5, 2, 0, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 1, 0, 3, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 1, 1, 4, 10, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 1, 5, 5, 10, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 1, 0, 0, 32, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 1, 1, 2, 32, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 1, 5, 4, 32, 10));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 0, 1, 0, -1, 5, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 1, 1, 1, -1, 33, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 1, 5, -1, 50, 10));
}

SEQAN_DEFINE_TEST(test_find_journaled_string_tree_fuzzy_occ_shift_and)
{
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, -1, 10, 10));;
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, 0, 11, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, 1, 12, 13));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, 2, 13, 9));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, 3, 14, 5));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, 4, 15, 21));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 5, 5, 16, 18));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, -1, 0, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, 0, 1, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, 1, 2, 13));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, 2, 3, 9));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, 3, 4, 5));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, 4, 5, 21));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 4, 6, 5, 6, 18));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, -1, 20, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, 0, 41, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, 1, 39, 13));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, 2, 85, 9));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, 3, 76, 5));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, 4, 50, 21));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 5, 5, 80, 18));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, -1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, 0, 41, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, 1, 2, 13));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, 2, 14, 9));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, 3, 18, 5));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, 4, 10, 21));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 3, 5, 6, 5, 34, 18));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, -1, 20, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, 0, 41, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, 1, 39, 13));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, 2, 85, 9));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, 3, 76, 5));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, 4, 50, 21));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 5, 5, 80, 18));

    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, -1, 0, 10));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, 0, 41, 11));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, 1, 2, 13));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, 2, 14, 9));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, 3, 18, 5));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, 4, 10, 21));
    SEQAN_ASSERT(_configureTest(char(), seqan::ShiftAnd(), 9, 5, 6, 5, 34, 18));
}

#endif  // EXTRAS_TESTS_FIND_JOURNALED_STRING_TREE_TEST_FIND_JOURNALED_STRING_TREE_SHIFT_AND_H_
