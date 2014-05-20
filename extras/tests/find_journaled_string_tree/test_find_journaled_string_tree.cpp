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
// Test suite for data parallel finder.
// ==========================================================================

//#define DEBUG_DATA_PARALLEL

#include <seqan/basic.h>

#include "../journaled_string_tree/test_journaled_string_tree_mock_generator.h"
#include "test_find_journaled_string_tree_compare_results.h"

#include "test_find_journaled_string_tree_simple.h"
#include "test_find_journaled_string_tree_horspool.h"
#include "test_find_journaled_string_tree_shift_and.h"
#include "test_find_journaled_string_tree_shift_or.h"
#include "test_find_journaled_string_tree_myers.h"

SEQAN_BEGIN_TESTSUITE(find_journaled_string_tree)
{
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ref_occ_simple);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_snp_occ_simple);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ins_occ_simple);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_del_occ_simple);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_fuzzy_occ_simple);

    SEQAN_CALL_TEST(test_find_journaled_string_tree_ref_occ_horspool);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_snp_occ_horspool);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ins_occ_horspool);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_del_occ_horspool);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_fuzzy_occ_horspool);

    SEQAN_CALL_TEST(test_find_journaled_string_tree_ref_occ_shift_and);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_snp_occ_shift_and);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ins_occ_shift_and);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_del_occ_shift_and);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_fuzzy_occ_shift_and);

    SEQAN_CALL_TEST(test_find_journaled_string_tree_ref_occ_shift_or);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_snp_occ_shift_or);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ins_occ_shift_or);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_del_occ_shift_or);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_fuzzy_occ_shift_or);

    SEQAN_CALL_TEST(test_find_journaled_string_tree_ref_occ_myers_ukkonen_infix_0_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ref_occ_myers_ukkonen_infix_3_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_snp_occ_myers_ukkonen_infix_0_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_snp_occ_myers_ukkonen_infix_3_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ins_occ_myers_ukkonen_infix_0_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_ins_occ_myers_ukkonen_infix_3_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_del_occ_myers_ukkonen_infix_0_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_del_occ_myers_ukkonen_infix_3_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_fuzzy_occ_myers_ukkonen_infix_0_errors);
    SEQAN_CALL_TEST(test_find_journaled_string_tree_fuzzy_occ_myers_ukkonen_infix_3_errors);

}
SEQAN_END_TESTSUITE
