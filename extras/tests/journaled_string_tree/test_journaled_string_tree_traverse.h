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
// Implements test for the find methods.
// ==========================================================================
#ifndef EXTRAS_TESTS_DATA_PARALLEL_TEST_DATA_PARALLEL_FIND_H_
#define EXTRAS_TESTS_DATA_PARALLEL_TEST_DATA_PARALLEL_FIND_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/journaled_string_tree.h>
#include <seqan/find_journaled_string_tree.h>

#include "test_journaled_string_tree_mock_generator.h"

// Global variable to read the config file.
seqan::DataParallelTestConfig<char> testConfig;

namespace seqan
{

struct TestDataParallel_;
typedef Tag<TestDataParallel_> TestDummy;

template <typename TNeedle>
class Pattern<TNeedle, TestDummy>
{
public:
    TNeedle _needle;

    // TODO(rmaerker): Test different pattern sizes.
    Pattern() : _needle("xxx")
    {}
};

template <typename TNeedle>
inline TNeedle & needle(Pattern<TNeedle, TestDummy> & pattern)
{
    return pattern._needle;
}

template <typename TNeedle>
inline TNeedle const & needle(Pattern<TNeedle, TestDummy> const & pattern)
{
    return pattern._needle;
}

template <typename TTraverser>
struct StringTreeTraversalTester
{
    TTraverser * _traverserPtr;

    StringTreeTraversalTester(TTraverser & other) : _traverserPtr(&other)
    {}
};

template <typename TTraverser, typename TTag>
inline typename ComputeState<TTraverser>::Type
compute(StringTreeTraversalTester<TTraverser> & /*tester*/,
        TTag const & /*tag*/)
{
    typedef typename ComputeState<TTraverser>::Type TCallState;
    return TCallState(true, 1);
}

template <typename TValue>
struct SequenceAppender
{

    StringSet<String<TValue> > _processedSeq;

    SequenceAppender(unsigned size)
    {
        resize(_processedSeq, size, Exact());
    }

    template <typename TTraverser>
    void operator()(TTraverser & traverser)
    {
        typedef typename Position<TTraverser>::Type TPosVec;
#ifdef TEST_DEBUG_OUTPUT
        std::cerr << "position(traverser): ";
#endif
        TPosVec posVec = position(traverser);
        for (unsigned i = 0; i < length(posVec); ++i)
        {
#ifdef TEST_DEBUG_OUTPUT
            std::cerr << "("<< posVec[i].i1 << ", " <<  posVec[i].i2 << ")" << "; ";
#endif
            appendValue(_processedSeq[posVec[i].i1],
                        value(journalData(container(traverser)), posVec[i].i1)[posVec[i].i2]);
        }
#ifdef TEST_DEBUG_OUTPUT
        std::cerr << std::endl;
#endif
    }
};

}

template <typename TTestSeq, typename TCompareSeq, typename TSize>
bool compareResults(TTestSeq const & testSeq, TCompareSeq const & compSeq, TSize const & windowSize)
{
    SEQAN_ASSERT_EQ(length(testSeq), length(compSeq));
    for (unsigned i = 0; i < length(testSeq); ++i)
        if (isNotEqual(testSeq[i], prefix(compSeq[i], length(compSeq[i]) - (windowSize -1))))
            return false;
    return true;
}


template <typename TMock, typename TTester, typename TSize>
void _printDebugInfo(TMock const & mockGen, TTester const & dpTester, TSize const & windowSize)
{
    std::cerr << "Host: " << host(journalData(mockGen._mock)) << std::endl;
    for (unsigned i = 0; i < length(dpTester._processedSeq); ++i)
        std::cerr << "Compare 1: " << dpTester._processedSeq[i] << "\nCompare 2: "
        << prefix(journalData(mockGen._mock)[i], length(mockGen._seqData[i]) - (windowSize - 1)) << "\n"
        << std::endl;
}

bool _runTestForConfiguration(unsigned posConf,
                              unsigned varConf,
                              unsigned covConf,
                              unsigned refLength,
                              unsigned windowSize,
                              seqan::StringTreeDefault const & /*stringTreeTag*/)
{
    using namespace seqan;

    typedef String<char, Journaled<Alloc<>, SortedArray > > TJournalString;
    typedef Host<TJournalString>::Type THost;
    typedef String<MockVariantData<char> > TVarData;
    typedef String<String<bool, Packed<> > > TCovData;
    typedef SequenceAppender<char> TSequenceAppender;
    typedef MockGenerator_<unsigned, char> TMockGenerator;
    typedef TMockGenerator::TStringTree TStringTree;
    typedef Traverser<TStringTree, TraverserSpec<> > TTraverser;
    typedef StringTreeTraversalTester<TTraverser> TDummyCaller;

    TVarData varData;
    TCovData covData;
    testConfig.getTestConfiguration(varData, covData, posConf, varConf, covConf);

    // Initialize the mock generator.
    TMockGenerator mockGen;
    // Generate the mock for the current configuration.
    mockGen.generate(varData, covData, refLength);

    TSequenceAppender seqAppender(length(journalData(mockGen._mock)));
    TTraverser traverser(mockGen._mock, windowSize);
    TDummyCaller dummyCaller(traverser);
    traverse(traverser, dummyCaller, seqAppender);

#ifdef TEST_DEBUG_OUTPUT
    _printDebugInfo(mockGen, seqAppender, windowSize);
#endif

    return compareResults(seqAppender._processedSeq, mockGen._seqData, windowSize);
}

// ----------------------------------------------------------------------------
// Test all SNPs.
// ----------------------------------------------------------------------------

// Test all at position 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 5, 101, 30, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 5, 101, 30, seqan::StringTreeDefault()));
}

// ----------------------------------------------------------------------------
// Test all deletions.
// ----------------------------------------------------------------------------

// Test all beginning at position 0, all dels, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_1_5_journaled_string_tree)
{
//    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// ----------------------------------------------------------------------------
// Test all Insertions.
// ----------------------------------------------------------------------------

// Test all at position 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_0_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_1_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_2_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_3_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 5, 101, 30, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_4_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 5, 101, 30, seqan::StringTreeDefault()));
}

// ----------------------------------------------------------------------------
// Test special variant comobinations.
// ----------------------------------------------------------------------------

// Test all deletions with different size.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_5_3_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_5_3_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_5_3_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_5_3_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_5_3_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_5_3_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 5, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_6_4_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(6, 4, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(6, 4, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_7_4_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(7, 4, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(7, 4, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_9_5_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_9_5_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_9_5_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_finder_config_9_5_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 5, 101, 10, seqan::StringTreeDefault()));
}


SEQAN_DEFINE_TEST(test_journaled_string_tree_merge_point_stack)
{
    using namespace seqan;

    typedef String<bool, Packed<> > TCoverage;
    typedef MergePointMap_<unsigned, TCoverage> TMergePointStore;


    TMergePointStore store;
    TCoverage test1;
    resize(test1, 10, false, Exact());
    test1[0] = true;
    test1[5] = true;
    TCoverage test2;
    resize(test2, 10, false, Exact());
    test2[0] = true;
    test2[6] = true;
    test2[8] = true;
    TCoverage test3;
    resize(test3, 10, false, Exact());
    test3[1] = true;
    test3[4] = true;
    test3[9] = true;

    push(store, 10, test1);
    push(store, 30, test2);
    push(store, 20, test3);
    push(store, 25, test2);

    TCoverage result;
    resize(result, 10, false, Exact());

    for (unsigned i = 0; i < length(store._mergePoints); ++i)
        std::cerr << "What: " << store._mergePoints[i] << std::endl;

    _syncToMergePoint(result, store, 5, MergePointSyncResize());
    std::cerr << result << std::endl;

    for (unsigned i = 0; i < length(store._mergePoints); ++i)
        std::cerr << "What: " << store._mergePoints[i] << std::endl;

    _syncToMergePoint(result, store, 15, MergePointSyncResize());
    std::cerr << result << std::endl;

    for (unsigned i = 0; i < length(store._mergePoints); ++i)
        std::cerr << "What: " << store._mergePoints[i] << std::endl;

    _syncToMergePoint(result, store, 25, MergePointSyncResize());
    std::cerr << result << std::endl;

    for (unsigned i = 0; i < length(store._mergePoints); ++i)
        std::cerr << "What: " << store._mergePoints[i] << std::endl;

    _syncToMergePoint(result, store, 60, MergePointSyncResize());
    std::cerr << result << std::endl;

    for (unsigned i = 0; i < length(store._mergePoints); ++i)

        std::cerr << "What: " << store._mergePoints[i] << std::endl;
}

#endif  // EXTRAS_TESTS_DATA_PARALLEL_TEST_DATA_PARALLEL_FIND_H_
