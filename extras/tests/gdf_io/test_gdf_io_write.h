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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef EXTRAS_TESTS_JSEQ_IO_TEST_JSEQ_IO_WRITE_H_
#define EXTRAS_TESTS_JSEQ_IO_TEST_JSEQ_IO_WRITE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/journaled_set.h>
#include <seqan/gdf_io.h>

SEQAN_DEFINE_TEST(test_jseq_io_write)
{
    using namespace seqan;

    typedef String<char, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournalString;
    typedef Host<TJournalString>::Type THost;
    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournalSet;

    THost ref = "A Man, A Plan, A Canal - Panama";

    CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/extras/tests/gdf_io/");
    CharString fileName = "goldGDF.gdf";
    unsigned refHash = 1234567890;


    TJournalString str1;
    setHost(str1, ref);
    erase(str1, 4, 6);
    insert(str1, 6, "n");
                        // 0123456789012345678901234567890
                        //"A Ma An Plan, A Canal - Panama";
    SEQAN_ASSERT_EQ(str1, "A Ma An Plan, A Canal - Panama");

    TJournalString str2;
    setHost(str2, ref);
    insert(str2, 8, "n");
    insert(str2, 17, "n");
    insert(str2, 20, "h");
    assignValue(str2, 29, 'o');
    assignValue(str2, 31, 'o');

    SEQAN_ASSERT_EQ(str2, "A Man, An Plan, An Chanal - Ponoma");

    TJournalString str3;
    setHost(str3, ref);
    erase(str3, 4, 6);
    insert(str3, 14, "n");
    insert(str3, 17, "h");
    assignValue(str3, 26, 'o');
    erase(str3, 28);
                        // 012345678901234567890123456789
                        //"A Ma A Plan, An Chanal - Panma";
    SEQAN_ASSERT_EQ(str3, "A Ma A Plan, An Chanal - Ponma");

    TJournalSet set;
    setHost(set, ref);
    appendValue(set, str1);
    appendValue(set, str2);
    appendValue(set, str3);

    JSeqHeader jseqHeader;
    jseqHeader._refInfos._refId = "Test Sequence";
    jseqHeader._refInfos._refFile = "File://Somewhere";
    jseqHeader._refInfos._refHash = refHash;

    jseqHeader._fileInfos._majorFileId = 2u;
    jseqHeader._fileInfos._minorFileId = 5u;
    jseqHeader._fileInfos._byteOrder = SystemsByteOrder::IS_LITTLE_ENDIAN();
    jseqHeader._fileInfos._blockSize = 1000000u;
    jseqHeader._fileInfos._snpCompression = false;

    appendValue(jseqHeader._nameStore, "seq1");
    appendValue(jseqHeader._nameStore, "seq2");
    appendValue(jseqHeader._nameStore, "seq3");

    appendValue(jseqHeader._headerRecords, JSeqHeaderRecord("This", "is"));
    appendValue(jseqHeader._headerRecords, JSeqHeaderRecord("a", "test"));

    CharString file = filePath;
    append(file, fileName);
    std::ofstream outStream;
    outStream.open(toCString(file), std::ios_base::out);
    if (!outStream.good())
    {
        SEQAN_ASSERT_FAIL("Cannot open file: %s", toCString(file));
    }

    write(outStream, set, jseqHeader, JSeq());
    outStream.close();
}

#endif  // EXTRAS_TESTS_JSEQ_IO_TEST_JSEQ_IO_WRITE_H_
