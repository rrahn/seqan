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
// Tests the partial sum manager.
// ==========================================================================
#ifndef SANDBOX_RMAERKER_TESTS_TEST_DELTA_INDEX_TEST_PARTIAL_SUM_MANAGER_H_
#define SANDBOX_RMAERKER_TESTS_TEST_DELTA_INDEX_TEST_PARTIAL_SUM_MANAGER_H_

#include <seqan/basic.h>
#include <seqan/index_journaled.h>

SEQAN_DEFINE_TEST(test_partial_sum_manager_insert_simple)
{
    using namespace seqan;

    PartialSumManager<int, Simple> psMgr;
    insert(psMgr, 10, 4);

    insert(psMgr, 3, -3);
    insert(psMgr, 18, -7);
    insert(psMgr, 3, 3);

    SEQAN_ASSERT_EQ(length(psMgr._partialSumTable), 3u);
    SEQAN_ASSERT_EQ(psMgr._partialSumTable[0].i1, 3);
    SEQAN_ASSERT_EQ(psMgr._partialSumTable[0].i2, 0);
    SEQAN_ASSERT_EQ(psMgr._partialSumTable[1].i1, 10);
    SEQAN_ASSERT_EQ(psMgr._partialSumTable[1].i2, 4);
    SEQAN_ASSERT_EQ(psMgr._partialSumTable[2].i1, 18);
    SEQAN_ASSERT_EQ(psMgr._partialSumTable[2].i2, -7);
}

SEQAN_DEFINE_TEST(test_partial_sum_manager_partial_sum_simple)
{
    using namespace seqan;

    PartialSumManager<int, Simple> psMgr;
    insert(psMgr, 10, 4);
    insert(psMgr, 3, -3);
    insert(psMgr, 18, -7);

    SEQAN_ASSERT_EQ(partialSum_(psMgr, 0), 0);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 1), 0);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 2), 0);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 3), -3);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 9), -3);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 10), 1);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 17), 1);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 18), -6);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 20), -6);

    insert(psMgr, 3, 3);

    SEQAN_ASSERT_EQ(partialSum_(psMgr, 0), 0);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 3), 0);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 9), 0);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 10), 4);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 17), 4);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 18), -3);
    SEQAN_ASSERT_EQ(partialSum_(psMgr, 20), -3);
}

#endif  // SANDBOX_RMAERKER_TESTS_TEST_DELTA_INDEX_TEST_PARTIAL_SUM_MANAGER_H_
