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
// Tests for reading jseq data.
// ==========================================================================

#ifndef EXTRAS_TESTS_JSEQ_IO_TEST_JSEQ_IO_READ_H_
#define EXTRAS_TESTS_JSEQ_IO_TEST_JSEQ_IO_READ_H_

#include <seqan/basic.h>

SEQAN_DEFINE_TEST(test_jseq_io_read)
{
    using namespace seqan;

    typedef String<char, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournalString;
    typedef Host<TJournalString>::Type THost;
    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournalSet;

    THost ref = "A Man, A Plan, A Canal - Panama";

    CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/extras/tests/gdf_io/");
    CharString fileName = "goldGDF.gdf";
//    unsigned refHash = 1234567890;

    CharString file = filePath;
    append(file, fileName);
    std::ifstream inputStream;
    inputStream.open(toCString(file), std::ios_base::in);
    if (!inputStream.good())
        SEQAN_ASSERT_FAIL("Cannot open file: %s", toCString(file));

    JSeqHeader jseqHeader;
    RecordReader<std::ifstream, SinglePass<> > reader(inputStream);
    TJournalSet journalSet;
    setHost(journalSet, ref);

    read(journalSet, jseqHeader, reader, JSeq());
    inputStream.close();

    // Test mandatory fields regarding the reference.
    SEQAN_ASSERT_EQ(jseqHeader._refInfos._refId , "Test Sequence");
    SEQAN_ASSERT_EQ(jseqHeader._refInfos._refFile , "File://Somewhere");
    SEQAN_ASSERT_EQ(jseqHeader._refInfos._refHash , 1234567890u);

    // Test mandatory fields regarding the file itself.
    SEQAN_ASSERT_EQ(jseqHeader._fileInfos._majorFileId , 0u);
    SEQAN_ASSERT_EQ(jseqHeader._fileInfos._minorFileId , 1u);
    SEQAN_ASSERT_EQ(jseqHeader._fileInfos._byteOrder , true);
    SEQAN_ASSERT_EQ(jseqHeader._fileInfos._blockSize , 1000000u);
    SEQAN_ASSERT_EQ(jseqHeader._fileInfos._snpCompression , false);

    // Test optional fields
    SEQAN_ASSERT_EQ(jseqHeader._headerRecords[0]._key, "This");
    SEQAN_ASSERT_EQ(jseqHeader._headerRecords[0]._value, "is");
    SEQAN_ASSERT_EQ(jseqHeader._headerRecords[1]._key, "a");
    SEQAN_ASSERT_EQ(jseqHeader._headerRecords[1]._value, "test");

    // Cannot write snps in compressed form if we use a char as alphabet.
    SEQAN_ASSERT_EQ(length(journalSet), 3u);
    SEQAN_ASSERT_EQ(host(journalSet[0]), ref);
    SEQAN_ASSERT_EQ(journalSet[0], "A Ma An Plan, A Canal - Panama");
    SEQAN_ASSERT_EQ(host(journalSet[1]), ref);
    SEQAN_ASSERT_EQ(journalSet[1], "A Man, An Plan, An Chanal - Ponoma");
    SEQAN_ASSERT_EQ(host(journalSet[2]), ref);
    SEQAN_ASSERT_EQ(journalSet[2], "A Ma A Plan, An Chanal - Ponma");
}

#endif  // EXTRAS_TESTS_JSEQ_IO_TEST_JSEQ_IO_READ_H_
