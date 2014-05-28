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
// Implements the version to compare results.
// ==========================================================================

#ifndef EXTRAS_TESTS_FIND_JOURNALED_STRING_TREE_TEST_FIND_JOURNALED_STRING_TREE_COMPARE_RESULTS_H_
#define EXTRAS_TESTS_FIND_JOURNALED_STRING_TREE_TEST_FIND_JOURNALED_STRING_TREE_COMPARE_RESULTS_H_

#include <algorithm>
#include <seqan/basic.h>
#include <seqan/journaled_set.h>
#include <seqan/journaled_string_tree.h>
#include <seqan/find_journaled_string_tree.h>
#include <seqan/random.h>

#include "../journaled_string_tree/test_journaled_string_tree_mock_generator.h"

// Global variable to read the config file.
seqan::DataParallelTestConfig<char>         testFinderConfigChar;
seqan::DataParallelTestConfig<seqan::Dna>   testFinderConfigDna;

seqan::Rng<seqan::MersenneTwister> rng(42);


template <typename TPosition>
struct FindResultsTester_
{
    typedef seqan::String<TPosition> THitString;

    seqan::StringSet<THitString> _hitStringSet;

    FindResultsTester_(unsigned size)
    {
        resize(_hitStringSet, size, seqan::Exact());
    }

    template <typename TTraverser>
    inline void operator()(TTraverser & traverser)
    {
        typedef typename seqan::Positions<TTraverser>::Type TPositionVec;

        TPositionVec posVec = positions(traverser);
        for (unsigned i = 0; i < length(posVec); ++i)
            appendValue(_hitStringSet[posVec[i].i1], posVec[i].i2);
    }
};

template <typename TMock, typename TJst, typename TNdl, typename TTheirHits, typename TMyHits, typename TIdx>
void
_printDebugInfo(TMock const & mockGen, TJst const & jst, TNdl const & ndl,
                TTheirHits const & theirHits, TMyHits const & myHits, TIdx const & idx)
{
    std::cout << "\nTheir Hits for " << idx<< std::endl;

    for (unsigned j = 0; j < length(theirHits); ++j)
        std::cout << theirHits[j] << ", ";

    std::cout << "\nMy Hits for " << idx << std::endl;
    for (unsigned j = 0; j < length(myHits[idx]); ++j)
        std::cout << myHits[idx][j] << ", ";

    std::cout << "\n\nSequences"<< std::endl;
    std::cout << "Pattern: " << ndl << std::endl;
    std::cout << "\nRef: " << host(mockGen._seqData) << "\n" << std::endl;

    for (unsigned i = 0; i < length(mockGen._seqData); ++i)
        std::cout << "Gen: " << mockGen._seqData[i] << std::endl;

    std::cout << "\n" << std::endl;
    for (unsigned i = 0; i < length(mockGen._seqData); ++i)
        std::cout << "JST: " << stringSet(jst)[i] << std::endl;
}

template <typename TMockGenerator, typename TAlphabet, typename TAlgorithm>
bool _runTest(TMockGenerator & mockGen,
              TAlphabet /*alphabet*/,
              TAlgorithm /*algo*/,
              int seqBase,
              unsigned patternPos,
              unsigned patternLength,
              int /*numErrors*/,
              unsigned blockSize)
{
    using namespace seqan;

    typedef String<TAlphabet> TNeedle;
//    typedef DeltaStore<unsigned, TAlphabet> TDeltaStore;
    typedef JournaledStringTree<DeltaMap<unsigned, TAlphabet>, StringTreeDefault> TContainer;

    typedef Pattern<TNeedle, TAlgorithm> TPattern;
    typedef Finder_<TContainer, TPattern, Jst<> > TFinder;
    typedef typename Position<TContainer>::Type TPosition;

    typedef typename TMockGenerator::TJournalSet TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
//    typedef typename Iterator<TJournalString, Standard>::Type TJournalIter;

    typedef ContainerView<TJournalString> TJournalView;
    typedef Finder<TJournalView> TStandardFinder;

    TNeedle ndl;
    if (seqBase == -1)
        ndl = infix(host(mockGen._seqData), patternPos, patternPos + patternLength);
    else
        ndl = infix(mockGen._seqData[seqBase], patternPos, patternPos + patternLength);

    TPattern pattern;
    setHost(pattern, ndl);

    FindResultsTester_<TPosition> delegate(length(mockGen._seqData));
    TContainer jst(host(mockGen._seqData), mockGen._varStore);

    if (blockSize > 0u)
        setBlockSize(jst, blockSize);

    TFinder finder(jst);
    find(finder, pattern, delegate);

    // Search over each sequence separately and find the pattern.
    for (unsigned i = 0; i < length(mockGen._seqData); ++i)
    {
        TJournalView jView(mockGen._seqData[i]);
        TStandardFinder stdFinder(jView);
        String<unsigned> compareHits;
        while (find(stdFinder, pattern))
        {
            if (IsSameType<TAlgorithm, Simple>::VALUE)
                appendValue(compareHits, position(stdFinder));
            else
                appendValue(compareHits, position(stdFinder) + (length(ndl) - 1));
        }
        if (length(compareHits) != length(delegate._hitStringSet[i]))
            return false;

        for (unsigned j = 0; j < length(compareHits); ++j)
            if (compareHits[j] != delegate._hitStringSet[i][j])
                return false;
    }
    return true;
}

// Test for approximate search
template <typename TMockGenerator, typename TAlphabet, typename TMyersSpec, typename THasState, typename TBeginSpec>
bool _runTest(TMockGenerator & mockGen,
              TAlphabet /*alphabet*/,
              seqan::Myers<TMyersSpec, THasState, TBeginSpec> /*algo*/,
              int seqBase,
              unsigned patternPos,
              unsigned patternLength,
              int numErrors,
              unsigned blockSize)
{
    using namespace seqan;

    typedef String<TAlphabet> TNeedle;
//    typedef DeltaStore<unsigned, TAlphabet> TDeltaStore;
    typedef JournaledStringTree<DeltaMap<unsigned, TAlphabet>, StringTreeDefault> TContainer;

    typedef Pattern<TNeedle, Myers<TMyersSpec, THasState, TBeginSpec> > TPattern;
    typedef Finder_<TContainer, TPattern, Jst<> > TFinder;
    typedef typename Position<TContainer>::Type TPosition;

    typedef typename TMockGenerator::TJournalSet TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
//    typedef typename Iterator<TJournalString, Standard>::Type TJournalIter;

    typedef ContainerView<TJournalString> TJournalView;
    typedef Finder<TJournalView> TStandardFinder;

    typedef Pdf<Uniform<unsigned> > TPdf;

    TNeedle ndl;
    if (seqBase == -1)
        ndl = infix(host(mockGen._seqData), patternPos, patternPos + patternLength);
    else
        ndl = infix(mockGen._seqData[seqBase], patternPos, patternPos + patternLength);

    // Alter the needle to search approximatively.
    for (int i = 0; i < std::abs(numErrors); ++i)
    {
        unsigned errorType = pickRandomNumber(rng, TPdf(0, 2));
        TAlphabet mut = static_cast<TAlphabet>(pickRandomNumber(rng, TPdf(65, 90)));
        unsigned mutPos = pickRandomNumber(rng, TPdf(0, length(ndl) - 1));
        switch(errorType)
        {
            case 0:
                ndl[mutPos] = mut;
                break;
            case 1:
                insertValue(ndl, mutPos, mut);
                break;
            case 2:
                erase(ndl, mutPos);
        }
    }

    TPattern pattern;
    setHost(pattern, ndl);

    FindResultsTester_<TPosition> delegate(length(mockGen._seqData));
    TContainer jst(host(mockGen._seqData), mockGen._varStore);
    if (blockSize > 0u)
        setBlockSize(jst, blockSize);
    TFinder finder(jst);
    find(finder, pattern, delegate, numErrors);

    // Search over each sequence separately and find the pattern.
    for (unsigned i = 0; i < length(mockGen._seqData); ++i)
    {
        TJournalView jView(mockGen._seqData[i]);

        TStandardFinder standardFinder(jView);

        String<unsigned> compareHits;
        while (find(standardFinder, pattern, numErrors))
            appendValue(compareHits, position(standardFinder));

        if (length(compareHits) != length(delegate._hitStringSet[i]))
        {
            _printDebugInfo(mockGen, jst, ndl, compareHits, delegate._hitStringSet, i);
            return false;
        }

        for (unsigned j = 0; j < length(compareHits); ++j)
            if (compareHits[j] != delegate._hitStringSet[i][j])
            {
                _printDebugInfo(mockGen, jst, ndl, compareHits, delegate._hitStringSet, i);
                return false;
            }
    }
    return true;
}

template <typename TAlgorithm>
bool _configureTest(seqan::Dna /*alphabet*/,
                    TAlgorithm /*algorithm*/,
                    unsigned posConf,
                    unsigned varConf,
                    unsigned covConf,
                    int seqBase,
                    unsigned patternPos,
                    unsigned patternLength,
                    int numErrors = 0,
                    unsigned blockSize = 0u)
{
    using namespace seqan;

    typedef String<MockVariantData<Dna> > TVarData;
    typedef String<String<bool, Packed<> > > TCovData;
    typedef MockGenerator_<unsigned, Dna> TMockGenerator;

    TVarData varData;
    TCovData covData;
    testFinderConfigDna.getTestConfiguration(varData, covData, posConf, varConf, covConf);

    // Initialize the mock generator.
    TMockGenerator mockGen;
    // Generate the mock for the current configuration.
    mockGen.generate(varData, covData, 110);

    return _runTest(mockGen, Dna(), TAlgorithm(), seqBase, patternPos, patternLength, numErrors, blockSize);
}

template <typename TAlgorithm>
bool _configureTest(char /*alphabet*/,
                    TAlgorithm /*algorithm*/,
                    unsigned posConf,
                    unsigned varConf,
                    unsigned covConf,
                    int seqBase,
                    unsigned patternPos,
                    unsigned patternLength,
                    int numErros = 0,
                    unsigned blockSize = 0u)
{
    using namespace seqan;

//    typedef String<char, Journaled<Alloc<>, SortedArray > > TJournalString;
//    typedef Host<TJournalString>::Type THost;
    typedef String<MockVariantData<char> > TVarData;
    typedef String<String<bool, Packed<> > > TCovData;

    typedef MockGenerator_<unsigned, char> TMockGenerator;

    TVarData varData;
    TCovData covData;
    testFinderConfigChar.getTestConfiguration(varData, covData, posConf, varConf, covConf);

    // Initialize the mock generator.
    TMockGenerator mockGen;
    // Generate the mock for the current configuration.
    mockGen.generate(varData, covData, 110);

    return _runTest(mockGen, char(), TAlgorithm(), seqBase, patternPos, patternLength, numErros, blockSize);
}

#endif  // EXTRAS_TESTS_FIND_JOURNALED_STRING_TREE_TEST_FIND_JOURNALED_STRING_TREE_COMPARE_RESULTS_H_
