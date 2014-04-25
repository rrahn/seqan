// ==========================================================================
//                               journaled_string_tree
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

#ifndef EXTRAS_TESTS_DATA_PARALLEL_TEST_DATA_PARALLEL_H_
#define EXTRAS_TESTS_DATA_PARALLEL_TEST_DATA_PARALLEL_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/journaled_string_tree.h>

SEQAN_DEFINE_TEST(test_journaled_string_tree_test)
{

    using namespace seqan;

    typedef String<Dna, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournalString;
    typedef typename Host<TJournalString>::Type THostString;
//    unsigned lastEntryId = 0;

                               //1         2
                     //012345678901234567890123456
    THostString ref = "ACGTACGTACGACACCGTACGTACGAC";
                     //ACGTACGTAC   T A  ACGTACGAC
                     //ACGTACGTAC   T A  ACGTACGAC
                     //0123456789   0 1  234567890
                     //---------xxxxx
                                  //xxxxx

// SNP at 23 -> find virPos 10!!!
    // find range: 20 - 28
    // for refPositions at this site



    TJournalString js;
    setHost(js, ref);

    erase(js, 15, 17);
    assignValue(js, 14, 'A');
    assignValue(js, 13, 'T');
    erase(js, 10, 13);

//    std::cerr << "JournalString: " << js << std::endl;
//    std::cerr << "hostToVirtualPos(9): " << hostToVirtualPosition(js, 9) << std::endl;
//    std::cerr << "hostToVirtualPos(10): " << hostToVirtualPosition(js, 10) << std::endl;
//    std::cerr << "hostToVirtualPos(11): " << hostToVirtualPosition(js, 11) << std::endl;
//    std::cerr << "hostToVirtualPos(12): " << hostToVirtualPosition(js, 12) << std::endl;
//    std::cerr << "hostToVirtualPos(13): " << hostToVirtualPosition(js, 13) << std::endl;
//    std::cerr << "hostToVirtualPos(14): " << hostToVirtualPosition(js, 14) << std::endl;
//    std::cerr << "hostToVirtualPos(15): " << hostToVirtualPosition(js, 15) << std::endl;
//    std::cerr << "hostToVirtualPos(16): " << hostToVirtualPosition(js, 16) << std::endl;
//    std::cerr << "hostToVirtualPos(17): " << hostToVirtualPosition(js, 17) << std::endl;
//    std::cerr << "hostToVirtualPos(18): " << hostToVirtualPosition(js, 18) << std::endl;

}

// A test for strings.
SEQAN_DEFINE_TEST(test_journaled_string_tree_util_extract_variant)
{
    using namespace seqan;

//    typedef String<Dna, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournalString;
//    typedef typename Iterator<TJournalString>::Type TJournalStingIterator;
//    typedef typename Host<TJournalString>::Type THostString;
//    unsigned lastEntryId = 0;
//    THostString ref = "ACGTACGTACGTACGTACGTACGT";
//
//    // Test insertion at beginning.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        insert(js, 0, "TTT");
//        insert(js, 10, "C");
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 0, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_INSERTION);
//        SEQAN_ASSERT_EQ(lastEntryId, 0u);
//        SEQAN_ASSERT_EQ(jsIt == begin(js, Standard()), true);
//        SEQAN_ASSERT_EQ(delSize, 0u);
//    }
//
//    // Test insertion in middle.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        insert(js, 10, "TTT");
//        insert(js, 20, "C");
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 10, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_INSERTION);
//        TJournalStingIterator compIt = begin(js, Standard()) + 10;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
//        SEQAN_ASSERT_EQ(delSize, 0u);
//    }
//
//    // Test insertion at end.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        insert(js, length(js), "TTT");
//        insert(js, 10, "C");
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, length(ref) +1, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_INSERTION);
//        TJournalStingIterator compIt = begin(js, Standard()) + length(ref) + 1;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
//        SEQAN_ASSERT_EQ(delSize, 0u);
//    }
//
//    // Test deletion at beginning.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        erase(js, 0, 4);
//        erase(js, 10, 12);
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 0, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_DELETION);
//        SEQAN_ASSERT_EQ(jsIt == begin(js, Standard()), true);
//        SEQAN_ASSERT_EQ(delSize, 4u);
//    }
//
//    // Test deletion in middle.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        erase(js, 10, 14);
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 10, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_DELETION);
//        TJournalStingIterator compIt = begin(js, Standard()) + 9;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentHostIt, *compIt._currentHostIt);
////        SEQAN_ASSERT_EQ(*jsIt._hostSegmentBegin, *compIt._hostSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._hostSegmentEnd), *(--compIt._hostSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 4u);
//
//        erase(js, 16, 17);
//        varType = _extractVariant(jsIt, delSize, js, 10, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_DELETION);
//        compIt = begin(js, Standard()) + 9;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentHostIt, *compIt._currentHostIt);
////        SEQAN_ASSERT_EQ(*jsIt._hostSegmentBegin, *compIt._hostSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._hostSegmentEnd), *(--compIt._hostSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 4u);
//
//        varType = _extractVariant(jsIt, delSize, js, 20, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_DELETION);
//        compIt = begin(js, Standard()) + 15;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentHostIt, *compIt._currentHostIt);
////        SEQAN_ASSERT_EQ(*jsIt._hostSegmentBegin, *compIt._hostSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._hostSegmentEnd), *(--compIt._hostSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 1u);
//    }
//
//    // Test deletion at end.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        erase(js, 10, 14);
//        erase(js, length(js) - 2, length(js));
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, length(ref) - 2, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_DELETION);
//        TJournalStingIterator compIt = begin(js, Standard()) + (length(js) - 1);
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentHostIt, *compIt._currentHostIt);
////        SEQAN_ASSERT_EQ(*jsIt._hostSegmentBegin, *compIt._hostSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._hostSegmentEnd), *(--compIt._hostSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 2u);
//    }
//
//    // Test SNP at beginning.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        assignValue(js, 0, 'T');
//        assignValue(js, 10, 'C');
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 0, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_SUBSTITUTION);
//        SEQAN_ASSERT_EQ(jsIt == begin(js, Standard()), true);
//        SEQAN_ASSERT_EQ(delSize, 1u);
//    }
//
//    // Test SNP in middle.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        assignValue(js, 0, 'T');
//        assignValue(js, 10, 'C');
//        assignValue(js, 15, 'G');
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 10, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_SUBSTITUTION);
//        TJournalStingIterator compIt = begin(js, Standard()) + 10;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentInsertionBufferIt, *compIt._currentInsertionBufferIt);
////        SEQAN_ASSERT_EQ(*jsIt._insertionBufferSegmentBegin, *compIt._insertionBufferSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._insertionBufferSegmentEnd), *(--compIt._insertionBufferSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 1u);
//    }
//
//    // Test SNP at End.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        assignValue(js, 0, 'T');
//        assignValue(js, 10, 'C');
//        assignValue(js, 15, 'G');
//        assignValue(js, 23, 'A');
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 23, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_SUBSTITUTION);
//        TJournalStingIterator compIt = begin(js, Standard()) + 23;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentInsertionBufferIt, *compIt._currentInsertionBufferIt);
////        SEQAN_ASSERT_EQ(*jsIt._insertionBufferSegmentBegin, *compIt._insertionBufferSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._insertionBufferSegmentEnd), *(--compIt._insertionBufferSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 1u);
//    }
//
//    // Test longer substitution at beginning.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        replace(js, 0, 3, "TCG");
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 0, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_SUBSTITUTION);
//        SEQAN_ASSERT_EQ(jsIt == begin(js, Standard()), true);
//        SEQAN_ASSERT_EQ(delSize, 3u);
//    }
//
//    // Test longer substitution in middle.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        replace(js, 0, 3, "TCG");
//        replace(js, 10, 12, "AA");
//        replace(js, 20, 21, "C");
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 10, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_SUBSTITUTION);
//        TJournalStingIterator compIt = begin(js, Standard()) + 10;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentInsertionBufferIt, *compIt._currentInsertionBufferIt);
////        SEQAN_ASSERT_EQ(*jsIt._insertionBufferSegmentBegin, *compIt._insertionBufferSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._insertionBufferSegmentEnd), *(--compIt._insertionBufferSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 2u);
//    }
//
//    // Test shadowed deletion at beginning.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        erase(js, 0, 3);
//        insert(js, 0, "AC");
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 0, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_INSERTION | ExtractDeltaType::VARIANT_DELETION);
//        SEQAN_ASSERT_EQ(jsIt == begin(js, Standard()), true);
//        SEQAN_ASSERT_EQ(delSize, 3u);
//
//        insert(js, 0, "ACT");
//
//        varType = _extractVariant(jsIt, delSize, js, 0, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_INSERTION | ExtractDeltaType::VARIANT_DELETION);
//        SEQAN_ASSERT_EQ(jsIt == begin(js, Standard()), true);
//        SEQAN_ASSERT_EQ(delSize, 3u);
//    }
//
//    // Test shadowed deletion at beginning.
//    {
//        TJournalString js;
//        setHost(js, ref);
//
//        erase(js, 0, 3);
//        insert(js, 0, "AC");  //-1
//
//        erase(js, 10, 15);
//        insert(js, 10, "CGT");
//        insert(js, 12, "AAA");
//        insert(js, 12, "TTT");
//        erase(js, 19, 21);
//        assignValue(js, 20, 'G');
//
//        unsigned delSize = 0;
//        TJournalStingIterator jsIt;
//        ExtractDeltaType::TValue varType = _extractVariant(jsIt, delSize, js, 11, lastEntryId);
//
//        SEQAN_ASSERT_EQ(varType, ExtractDeltaType::VARIANT_INSERTION | ExtractDeltaType::VARIANT_DELETION);
//        TJournalStingIterator compIt = begin(js, Standard()) + 10;
//        SEQAN_ASSERT_EQ(*jsIt._journalEntriesIterator, *compIt._journalEntriesIterator);
////        SEQAN_ASSERT_EQ(*jsIt._currentInsertionBufferIt, *compIt._currentInsertionBufferIt);
////        SEQAN_ASSERT_EQ(*jsIt._insertionBufferSegmentBegin, *compIt._insertionBufferSegmentBegin);
////        SEQAN_ASSERT_EQ(*(--jsIt._insertionBufferSegmentEnd), *(--compIt._insertionBufferSegmentEnd));
//        SEQAN_ASSERT_EQ(delSize, 7u);
//    }

}

#endif  // EXTRAS_TESTS_DATA_PARALLEL_TEST_DATA_PARALLEL_H_
