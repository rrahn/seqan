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
// Unit tests for the journaled index implementation.
// ==========================================================================

#ifndef SANDBOX_RMAERKER_TESTS_INDEX_JOURNALED_TEST_INDEX_JOURNALED_BASE_H_
#define SANDBOX_RMAERKER_TESTS_INDEX_JOURNALED_TEST_INDEX_JOURNALED_BASE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index_journaled.h>

SEQAN_DEFINE_TEST(test_delta_index_impl_construction)
{
    using namespace seqan;

//    DnaString text = "ACGTACGAGCGGGACTGACGAGCGAGTTTTATCGGATCGGACTAGCGATCGACGTAGCGGAGCTACGGACTACGTAGAGCTGACGTACGGACGACGGACTGACGGACTAGCGACGTAGCGACGTTTAGCGAGAAAAATCGAGCGAAAATCGAGCTACGAGCTACCGGAGCTACGGACTACGTAGCGGAGCTACGGACTACGTAGCGGAGCTACGGACTACGTAGCGGAGCTACGGACT";
//    DnaString text = "TGTGGACGGCCAGTGGTGGGTTGCTACACCCCTGCGGCAACGTTGAAGCTCCTGGATTACACTGGCTGGATCTAAGCCGTGACACCCGTCATACTCCATAACCGTCTGTAACTCACGGCTTGTTCTGGACTGGATTGCCATTCTCTCAGAGTATTATGCAGGCCGGCGTACGGGTCCCATATAAACCTGTCATAGCTTACCTGACTCTACTTGGAAATGTGGCTAGGTCTTTGCCCACGCACCTAATCGGTCCTCGTTTGCTTTTTAGGACCCGATGAACTACAGAACACTGCAAGAATCTCTACCTGCTTTACAAAGTGCTGGATCCTATTCCAGCGGGATCTTTTATCTAAACACGATGAGAGGAGTATTCGTCAGGCCACATAGCTTTCTTGTTCTGATCGGAACGATCGTTGGCGCCCGACCCCCCGATTCCATAGTGAGTTCTTCGTCCGAGCCATTGTATGCGAGATCGATAGACTGATAGGGGATGCAGTATATCCCTGGATACAATAGACGCACAGGTTGGAATCCTAAGTGAAGTCGCGCGTCCGAACCCAGCTCTATTTTAGAGGTCATGGGTTCTGGTGCCCGCGAGCCGCGGAACCGATTAGGGGCATGTACAACAATATTTATTAGTCATCTTTCAGACACAATCTCCCAGCTCACTGGTATATAGTTCCTGCTATAATTAGCCTCCCTCATAAGTTGCACTACTTCAGCGTCCCAAATGCACCCTTACCACGAAGACAGGATTGTCCGATCCCATATTACGACCTTGGCAGGGGGTTCGCAAGTCCCACCCCAAACGATGCTGAAGGCTCAGGTTTCACAGGGACAAAAGCTTTAAACGCGAGTTCCCGCTCATAACCTGGACCGAATGCAGAATCATGCATCGTTCCACTGTGTTCGTGTCATCTAGGACGGGCGCAAAGGATATATAATTCAATTTTGAATACCTTATATTATTGTACACCTACCGGTCACCAGCCAACAATGT";
    DnaString text = "ACCGTACGAATCGA";

    String<Dna, Journaled<Alloc<>, SortedArray> > journalString(text);

    Index<DnaString, IndexEsa<> > esaIndex(text);

    Index<DnaString, IndexJournaled<IndexEsa<> > > deltaIndex(esaIndex, journalString);

    // We need dependencies between the construction scenarios.
    // At the moment the indices are set independently.
    indexCreate(deltaIndex, EsaSA());
    indexCreate(deltaIndex, EsaLcp());
    indexCreate(deltaIndex, EsaISA());

    std::cout << "\tPos\tSA\tiSA\tLCP\tText\n";
    for (unsigned i = 0; i < length(text); ++i)
    {
        std::cout << "\t" << i << "\t" << saAt(i, deltaIndex) << "\t" <<
                     invSaAt(i, deltaIndex) << "\t" << lcpAt(i, deltaIndex) << "\t" <<
                     suffix(journalString, saAt(i, deltaIndex)) << std::endl;
//        std::cout << "\t" << i << "\t" << saAt(i, deltaIndex) << "\t" <<
//                    invSaAt(i, deltaIndex) << "\t" << lcpAt(i, deltaIndex) << "\t" <<
//                    suffix(text, saAt(i, deltaIndex)) << std::endl;
    }

    _printDebug(indexText(deltaIndex));
    _printDebug(indexSA(deltaIndex));
    _printDebug(indexLcp(deltaIndex));

    insert(journalString, 7, "TAT");
    std::cout <<  journalString << std::endl;

//    std::cout << "\tPos\tSA\tLCP\tText\n";
//    for (unsigned i = 0; i < length(text); ++i)
//    {
//        std::cout << "\t" << i << "\t" << saAt(i, deltaIndex) << "\t" <<
//                    lcpAt(i, deltaIndex) << "\t" <<
//                    suffix(journalString, saAt(i, deltaIndex)) << std::endl;
////        std::cout << "\t" << i << "\t" << saAt(i, deltaIndex) << "\t" <<
////                    invSaAt(i, deltaIndex) << "\t" << lcpAt(i, deltaIndex) << "\t" <<
////                    suffix(text, saAt(i, deltaIndex)) << std::endl;
//    }

//    typedef typename Fibre<Index<DnaString, IndexJournaled<IndexEsa<> > >, FibreSA>::Type TSAFibre;
//    typedef typename Iterator<TSAFibre, Standard>::Type TSAIter;
//
//    TSAIter it = begin(indexSA(deltaIndex), Standard());
//    TSAIter itEnd = end(indexSA(deltaIndex), Standard());
//
////    std::cout << "The iter value: " << *it << std::endl;
//
//    insert(deltaIndex.psMgr, 5, 2);
//    insert(deltaIndex.psMgr, 10, -2);
////    std::cout << "The iter value: " << *it << std::endl;
//
//    while(it!=itEnd)
//    {
//        std::cout << "The iter value: " << *it << std::endl;
//        ++it;
//    }

    synchronize(deltaIndex);
    _printDebug(indexLcp(deltaIndex));
    std::cerr << "length(indexSA(deltaIndex)):  " << length(indexSA(deltaIndex)) << " length(indexLcp(deltaIndex)): " << length(indexLcp(deltaIndex)) << std::endl;
    std::cout << "\tPos\tSA\tLCP\tText\n";
    for (unsigned i = 0; i < length(journalString); ++i)
    {
        std::cout << "\t" << i << "\t" << saAt(i, deltaIndex /*indexSA(deltaIndex), i*/) << "\t" <<
                    lcpAt(i, deltaIndex) << "\t" <<
                    suffix(journalString, saAt(i, deltaIndex))  /*getValue(indexSA(deltaIndex), i))*/ << std::endl;
//        std::cout << "\t" << i << "\t" << saAt(i, deltaIndex) << "\t" <<
//                    invSaAt(i, deltaIndex) << "\t" << lcpAt(i, deltaIndex) << "\t" <<
//                    suffix(text, saAt(i, deltaIndex)) << std::endl;
    }
    _printDebug(indexLcp(deltaIndex));

    // Write your test here!

//    insert(deltaIndex, 10, "AA");

//    for (unsigned i = 0; i < length(text); ++i)
//    {
//        std::cout << saAt(i, value(deltaIndex.ref.ref)) << "\t";
//    }
//    std::cout << std::endl;
//    for (unsigned i = 0; i < length(text); ++i)
//    {
//        std::cout << lcpAt(i, value(deltaIndex.ref.ref)) << "\t";
//    }
//    std::cout << std::endl;
//
//    std::cout << "The inverse SA:\n";
//    for (unsigned i = 0; i < length(text); ++i)
//    {
//        std::cout << i << "\t";
//    }
//    std::cout << std::endl;
//    for (unsigned i = 0; i < length(text); ++i)
//    {
//        std::cout << deltaIndex.ref.invSa[i] << "\t";
//    }
//    std::cout << std::endl;

}

#endif  // SANDBOX_RMAERKER_TESTS_INDEX_JOURNALED_TEST_INDEX_JOURNALED_BASE_H_
