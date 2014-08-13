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
// Implements the general find method to search within the journal sequence.
// ==========================================================================

#ifndef EXTRAS_TOOLS_JSEQ_TOOLS_JSEQ_TOOLS_FIND_H_
#define EXTRAS_TOOLS_JSEQ_TOOLS_JSEQ_TOOLS_FIND_H_

#include "jst_bench_facade_header.h"

#include <seqan/journaled_string_tree.h>
#include <seqan/find_journaled_string_tree.h>

#include <seqan/random.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct HitCollector
{
    StringSet<String<unsigned> > hitPositionSet;
    String<int> _vpBlock;

    HitCollector(unsigned size)
    {
        resize(hitPositionSet, size, Exact());
        resize(_vpBlock, size, 0, Exact());
    }

    template <typename TTraverser>
    inline void operator()(TTraverser & finder)
    {
        typedef typename seqan::Positions<TTraverser>::Type TPosVec;

        SEQAN_OMP_PRAGMA(critical(insert))
        {
        TPosVec posVec = positions(finder);
        for (unsigned i = 0; i < length(posVec); ++i)
        {
            appendValue(hitPositionSet[posVec[i].i1], posVec[i].i2 + _vpBlock[posVec[i].i1]);
//            std::cout << "Hit: " << posVec[i] << std::endl;
        }
        }
    }

//    void print()
//    {
//        std::sort(begin(hitPosition), end(hitPosition));
//        for (unsigned i = 0; i < length(hitPosition); ++i)
//            std::cout << "Hit in seq: " << hitPosition[i].i1 << " at " << hitPosition[i].i2 << std::endl;
//    }
    void reinit()
    {
        for (unsigned i = 0; i < length(hitPositionSet); ++i)
            clear(hitPositionSet[i]);
    }

    template <typename TStream, typename TPatternId, typename TNameStore>
    void writeHits(TStream & stream, TPatternId const & patternId, TNameStore const & nameStore)
    {
        stream << "Pattern: " << patternId << '\n';
        for (unsigned i = 0; i < length(hitPositionSet); ++i)
        {
            stream << nameStore[i];
            for (unsigned j = 0; j < length(hitPositionSet[i]); ++j)
                stream << '\t' << hitPositionSet[i][j];
            stream << '\n';
        }
        stream << '\n';
    }
};

struct ParallelHitCollector_
{
    String<HitCollector> _hitCollectorSet;

    ParallelHitCollector_(unsigned numSeq, unsigned numThreads)
    {
        HitCollector tmp(numSeq);

        for (unsigned i = 0; i < numThreads; ++i)
            appendValue(_hitCollectorSet, tmp);
    }

    template <typename TExecutor>
    inline void operator()(TExecutor & executor)
    {
        // Delegate to the single threaded collector per thread.
        _hitCollectorSet[omp_get_thread_num()](executor);
    }
};


class ResultWriter_
{
public:

    ResultWriter_(CharString const & fileName)
    {
        oStream.open(toCString(fileName), std::ios_base::out);
    }

    template <typename TPatternId, typename TPatternSeq, typename THits>
    bool writeFormatted(TPatternId const & patternId,
                        TPatternSeq const & patternSeq,
                        THits const & patternHits)
    {
        if (!oStream.good())
            return false;

        oStream << "Pattern " << patternId << "\n" << patternSeq << std::endl;
        oStream << "BEGIN HITS\n";
        for (unsigned i = 0; i < length(patternHits); ++i)
        {
            oStream << "Hits for seq: " << i << ":\t";
            for (unsigned j = 0; j < length(patternHits[i]); ++j)
                oStream << patternHits[i][j] << ", ";
            oStream << "\n";
        }
        oStream << "\nEND HITS\n-----------------------------------------";
        return true;
    }

    ~ResultWriter_()
    {
        oStream.close();
    }

private:

    std::ofstream oStream;
};

// ----------------------------------------------------------------------------
// Struct FindContext_
// ----------------------------------------------------------------------------

struct FindContext_
{
    typedef Pair<int, unsigned> Hit;

    String<Hit> hitString;

    FindContext_() : hitString()
    {}

    void operator()(unsigned seq, unsigned pos)
    {
        appendValue(hitString, Pair<unsigned, unsigned>(seq,pos));
    }

    void print()
    {
        for (unsigned i = 0; i < length(hitString); ++i)
            std::cout << "Hit in seq: " << hitString[i].i1 << " at " << hitString[i].i2 << std::endl;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TResult, typename TFinder, typename TPattern, typename TSize, typename TSpec, typename THasState, typename TBeginSpec>
inline void
_findPatternForComparison(TResult & res,
                          TFinder &testFinder,TPattern & pattern,
                          TSize const & /*windowLength*/,
                          FindOptions const & findOptions,
                          Myers<TSpec, THasState, TBeginSpec> const & /*method*/)
{
    while(find(testFinder, pattern, findOptions.k))
        appendValue(res, position(testFinder));
}


template <typename TResult, typename TFinder, typename TPattern, typename TSize, typename TMethod>
inline void
_findPatternForComparison(TResult & res,
                          TFinder &testFinder,
                          TPattern & pattern,
                          TSize const & windowLength,
                          FindOptions const & /*findOptions*/,
                          TMethod const & /*method*/)
{
    while(find(testFinder, pattern))
    {
//                    std::cerr << " The position of the view: "  << position(testFinder) << std::endl;
        if (IsSameType<TMethod, Simple>::VALUE)
            appendValue(res, position(testFinder));
        else
            appendValue(res, position(testFinder) + (windowLength - 1));

    }
}

template <typename T>
void _printWithoutN(T const & seq)
{
    for (unsigned i = 0; i < length(seq); ++i)
        if (seq[i] != 'N')
            std::cout << seq[i];

    std::cout << std::endl;
}

template <typename TStringTree, typename TMethod>
void
_checkCorrectAlgorithm(TStringTree & stringTree,
                       FindOptions const & findOptions,
                       TMethod const &)
{
    typedef Pdf<Uniform<unsigned> > TPdf;
    typedef typename GetStringSet<TStringTree>::Type TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
    typedef typename Iterator<TJournalString>::Type TJournalIter;
    typedef typename Host<TJournalString>::Type THost;
    typedef typename Value<THost>::Type TAlphabet;

    typedef Pattern<THost, TMethod> TPattern;
    typedef Finder_<TStringTree, TPattern, Jst<TMethod> > TFinder;

    typedef ContainerView<TJournalString> TJournalStringView;
    typedef Finder<TJournalStringView> TTestFinder;

    Rng<MersenneTwister> rng(findOptions.seed);
    TPdf wlPdf(findOptions.wlRange.i1, findOptions.wlRange.i2);
//    TPdf deltaPdf(0, length(container(stringTree)) - 1);
    TPdf errorPdf(0, std::abs(findOptions.k));

    ResultWriter_ writer(findOptions.outputFile);

    double timeJournalGeneration = sysTime();

    TJournalSet compareSet;
    setHost(compareSet, host(stringSet(stringTree)));
    if (findOptions.numThreads > 1)
        adaptTo(compareSet, container(stringTree), 0, length(container(stringTree)), Parallel());
    else
        adaptTo(compareSet, container(stringTree), 0, length(container(stringTree)), Serial());

    unsigned errCount = 0;

//    _printWithoutN(host(compareSet));
//    if (findOptions.numThreads > 1)
//        requireJournal(stringTree, Parallel());
//    else
//        requireJournal(stringTree);
    if (findOptions.verbosity > 1)
        std::cout << "Time for Journal Generation: " << sysTime() - timeJournalGeneration << std::endl;

    unsigned numExp = _min(length(container(stringTree)), findOptions.patternCount);
    for (unsigned i = 1; i < numExp; ++i)
    {
        std::cerr << "\n########### Experiment " << i <<  " #############" << std::endl;
        std::cerr << "Seed: " << findOptions.seed << std::endl;
        unsigned windowLength = pickRandomNumber(rng, wlPdf);
        std::cerr << "WindowSize: " << windowLength << std::endl;

        unsigned deltaPos = i; //pickRandomNumber(rng, deltaPdf);
        std::cerr << "Delta: " << deltaPos << " -> " << container(stringTree)._keys[deltaPos] << std::endl;

        unsigned err = pickRandomNumber(rng, errorPdf);
        std::cerr << "Errors in needle: " << err << std::endl;

        unsigned relBeginPos = pickRandomNumber(rng, TPdf(0, windowLength + err -1));
        std::cerr << "RelBeginPos: " << relBeginPos << std::endl;

        if (testAllZeros(deltaCoverage(iter(container(stringTree), deltaPos, Standard()))))
        {
            std::cerr << "Zero coverage! Continue!" << std::endl;
            continue;
        }
        unsigned proxyId = bitScanForward(deltaCoverage(iter(container(stringTree), deltaPos, Standard())));
        std::cerr << "Proxy ID: " << proxyId << std::endl;
        TJournalString & js = value(compareSet, proxyId);
        TJournalIter itJs = begin(js, Standard());
        unsigned hostPos = container(stringTree)._keys[deltaPos];
        _mapHostToVirtual(itJs, js, container(stringTree), proxyId, hostPos);

        int beginPos = position(itJs) - relBeginPos;
        std::cerr << "VirtBeginPos: " << beginPos << std::endl;

        THost needle = infix(js, beginPos, beginPos + windowLength);
        std::cerr << "Needle: \t" << needle << std::endl;

        std::cerr << "Host-Region: \t" << infix(host(compareSet), hostPos - relBeginPos, hostPos - relBeginPos + windowLength) << std::endl;

        // Now we can also insert some errors into the needle
        for (unsigned j = 0; j < err; ++j)
        {
            // First select the type of
            unsigned errType = pickRandomNumber(rng, TPdf(0, 2));
            unsigned errPos = pickRandomNumber(rng, TPdf(0, length(needle) - 1));
            TAlphabet repChar;
            switch(errType)
            {
                case 0:
                    repChar = static_cast<TAlphabet>(pickRandomNumber(rng, TPdf(0, ValueSize<TAlphabet>::VALUE)));
                    needle[errPos] = repChar;
                    break;
                case 1:
                    repChar = static_cast<TAlphabet>(pickRandomNumber(rng, TPdf(0, ValueSize<TAlphabet>::VALUE)));
                    insertValue(needle, errPos, repChar);
                    break;
                case 2:
                    erase(needle, errPos);
            }
        }

        TPattern pattern;
        setHost(pattern, needle);

        reinit(stringTree);
        TFinder finder(stringTree);
        ParallelHitCollector_ delegateParallel(length(stringSet(stringTree)), findOptions.numThreads);

        if (findOptions.numThreads > 1)
            find(finder, pattern, delegateParallel, findOptions.k, Parallel());
        else
            find(finder, pattern, delegateParallel, findOptions.k, Serial());

        std::stringstream patternLable;
        patternLable << "Needle Information\t";
        patternLable << "WindowSize: " << windowLength << "\t";
        patternLable << "Delta: " << deltaPos << " -> " << container(stringTree)._keys[deltaPos] << "\t";
        patternLable << "BeginPos: " << beginPos << "\t";
        patternLable << "Proxy ID: " << proxyId << "\t";
        patternLable << "Needle: " << needle << "\t";
        patternLable << "Errors in needle: " << err << "\t";
        patternLable << "Errors allowed: " << findOptions.k << "\t";
        patternLable << "Seed: " << findOptions.seed << "\t";

        StringSet<String<unsigned> > dpPositions;
        resize(dpPositions, length(stringSet(stringTree)));

        for (unsigned q = 0; q < length(dpPositions); ++q)
        {
            for (unsigned p = 0; p < length(delegateParallel._hitCollectorSet); ++p)
                append(dpPositions[q], delegateParallel._hitCollectorSet[p].hitPositionSet[q]);
            if (!empty(dpPositions[q]))
                std::sort(begin(dpPositions[q], Standard()), end(dpPositions[q], Standard()));
        }

        if (!writer.writeFormatted(patternLable.str(), needle, dpPositions))
        {
            std::cerr << "Can't write to file: " << findOptions.outputFile << std::endl;
            return;
        }


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Search in each sequence separately for comparison.

        Pattern<THost, Horspool> exactPattern;
        setHost(exactPattern, needle);

        // Now we search within the journal sequences.
        StringSet<String<unsigned> > hits;
        resize(hits, length(compareSet));

        omp_set_num_threads(8);
        double counter = 0.0;
        Splitter<unsigned> jSetSplitter(0, length(hits), Parallel());
        SEQAN_OMP_PRAGMA(parallel for firstprivate(counter))
        for (unsigned jobId = 0; jobId < length(jSetSplitter); ++jobId)
        {
            counter = 0.0;
            unsigned currentBlock = 0;
            double intervalBorder = (double) (jSetSplitter[jobId + 1] - jSetSplitter[jobId]) * 0.1;
            printf("Search Sequential: 0%%");

            for (unsigned k = jSetSplitter[jobId]; k < jSetSplitter[jobId + 1]; ++k, ++counter)
            {
                // Search entire sequence.
                TJournalStringView jsView(compareSet[k]);
                TTestFinder testFinder(jsView);

                _findPatternForComparison(hits[k], testFinder, pattern, windowLength, findOptions, TMethod());

                if (counter >= intervalBorder)
                {
                    counter -= intervalBorder;
                    currentBlock += 10;
                    SEQAN_OMP_PRAGMA(critical(cout))
                    {
                        std::cout << " " << currentBlock << "%" << std::flush;
                    }
                }
            }
        }
        std::cout << " 100%" << std::endl;

        // Compare the results found by both methods.
        // We might found multiple hits?
        std::cerr << "--------------------- Check hits ---------------------" << std::endl;
        bool allGood = true;
        for (unsigned l = 0; l < length(hits); ++l)
        {
            // What can we do to check the correct result.
            String<Pair<int, int> > check;
            unsigned maxLength = _min(length(hits[l]), length(dpPositions[l]));

            for (unsigned m = 0; m < maxLength; ++m)
            {
                if (hits[l][m] != dpPositions[l][m])
                    appendValue(check, Pair<int, int>(hits[l][m], dpPositions[l][m]));
            }
            if (length(hits[l]) > maxLength)
                for (unsigned m = maxLength; m < length(hits[l]); ++m)
                    appendValue(check, Pair<int, int>(hits[l][m], -1));
            if (length(dpPositions[l]) > maxLength)
                for (unsigned m = maxLength; m < length(dpPositions[l]); ++m)
                    appendValue(check, Pair<int, int>(-1, dpPositions[l][m]));

            if (!empty(check))
            {
                std::cerr << "ERROR! Sequence " << l << ": ";
                for (unsigned m = 0; m < length(check); ++m)
                    std::cerr << " " << check[m];
                std::cerr << std::endl;
                allGood = false;
            }
        }
        if (allGood)
            std::cerr << "OK" << std::endl;
        else
            ++errCount;
    }

    std::cerr << "=================== Results ====================" << std::endl;
    std::cerr << "Experiments failed: \t\t" << errCount << " (" << numExp << ")" << std::endl;
    std::cerr << "Success rate: \t\t\t"     << static_cast<unsigned>((double) (numExp - errCount) / (double) numExp * 100) << " %"<< std::endl;
}


// ----------------------------------------------------------------------------
// Function _doFindPattern
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TNeedle, typename TSpec>
inline void _doFindPattern(FindContext_ & findContext,
                           StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournaledSet> > & journalSet,
                           Pattern<TNeedle, TSpec> & pattern,
                           FindOptions const & findOptions)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Position<TJournalString>::Type TPosition;
    typedef typename JournalType<TJournalString>::Type TJournalEntries;
    typedef typename GetCargoString_<TJournalEntries>::Type TEntryString;
    typedef typename Value<TEntryString>::Type TEntry;
    typedef typename Iterator<TEntryString, Standard>::Type TEntryIterator;
    typedef StringSet<TJournalString, Owner<JournaledSet> > TStringSet;
    typedef typename Host<TStringSet>::Type THost;
    typedef typename Iterator<THost, Standard>::Type THostIterator;
    typedef ContainerView<TJournalString> TJournalView;
    typedef ContainerView<THost> THostView;
    typedef typename Iterator<TStringSet, Standard>::Type TStringSetIterator;
//    typedef Segment<TJournalString, InfixSegment> TJournalInfix;


    // Start with the serial version here.

    // Variant one: We parse over each sequence of the journaled string and search for the pattern.
    double timeStart = sysTime();

//    // Parallelize with OpenMp
//    omp_set_num_threads(findOptions.numThreads);
//    Splitter<TStringSetIterator> setSplitter(begin(journalSet, Standard()), end(journalSet, Standard()), Parallel());
//
//    SEQAN_OMP_PRAGMA(parallel for)
//    for (unsigned job = 0; job < length(setSplitter); ++job)
//    {
//        Finder<TJournalString> journalFinder;
//        for (TStringSetIterator it = setSplitter[job]; it != setSplitter[job+1]; ++it)
//        {
//            setContainer(journalFinder, *it);
//            while(find(journalFinder, pattern))
//            {
//                SEQAN_OMP_PRAGMA(critical)
//                findContext(it - begin(journalSet, Standard()), position(journalFinder));
//            }
//            clear(journalFinder);
//        }
//    }
//    std::cout << "Time: " << sysTime() - timeStart << " s.\n";


    // Variant two: We search the reference sequence and then all diffs.
    unsigned overlapSize = length(host(pattern)) -1;

    // Parallelize with OpenMp
    omp_set_num_threads(findOptions.numThreads);

    Splitter<THostIterator> refSplitter(begin(host(journalSet), Standard()), end(host(journalSet), Standard()), Parallel());
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned job = 0; job < length(refSplitter); ++job)
    {
        printf("Thread %ld of %ld threads.\n", job, length(refSplitter));
        THostView tmpHostView;
        if (job == 0)
        {
            tmpHostView._begin = refSplitter[job];
            tmpHostView._end = refSplitter[job+1];
        }
        else
        {
            tmpHostView._begin = refSplitter[job] - overlapSize;
            tmpHostView._end = refSplitter[job+1];
        }

        Finder<THostView> refFinder(tmpHostView);
        while(find(refFinder, pattern))
        {
            SEQAN_OMP_PRAGMA(critical)
            findContext(-1, position(refFinder) + position(tmpHostView._begin, host(journalSet)));
        }
    }

    std::cout << "Break!\n";
    unsigned refHits = length(findContext.hitString);
    Splitter<TStringSetIterator> jobSplitter(begin(journalSet, Standard()), end(journalSet, Standard()), Parallel());
    SEQAN_OMP_PRAGMA(parallel for)
    for(unsigned job = 0u; job < length(jobSplitter); ++job)
    {
//        printf("Job %ld of %ld jobs.\n", job, length(jobSplitter));

        Finder<TJournalView> viewFinder;
        for (TStringSetIterator it = jobSplitter[job]; it != jobSplitter[job+1]; ++it)
        {

            // How can we check if a sequence covers a certain region in the reference.
            for (unsigned i = 0; i < refHits; ++i)
            {

                TPosition virtPos = hostToVirtualPosition(*it, findContext.hitString[i].i2);
                TEntry entry = *(findInJournalEntries(_journalEntries(*it), virtPos));
                if (entry.segmentSource == SOURCE_ORIGINAL)
                {
                    if (entry.virtualPosition + entry.length >= length(host(pattern)))
                    {
                        SEQAN_OMP_PRAGMA(critical)
                        findContext(it - begin(journalSet, Standard()), virtPos);
                    }
                }
            }
    //        std::cerr << "Processing sequence: " << it - begin(journalSet, Standard()) << std::endl;
            TEntryIterator entryIt = begin(_journalEntries(*it)._journalNodes, Standard());
            TJournalView tmpView(begin(*it), begin(*it));
            setContainer(viewFinder, tmpView);

            // Search the initial case.
    //        unsigned lastPhysicalPosOrigin = 0;
            if (entryIt->segmentSource == SOURCE_PATCH)
            {
                tmpView._end += entryIt->length + overlapSize;
                while(find(viewFinder, pattern))
                {
                    SEQAN_OMP_PRAGMA(critical)
                    findContext(it - begin(journalSet, Standard()), position(viewFinder) + position(haystack(viewFinder)._begin));
                }
            }
    //        else
    //            lastPhysicalPosOrigin += entryIt->physicalPosition + entryIt->length;

            haystack(viewFinder)._begin += (entryIt->length - overlapSize);
            hostIterator(viewFinder).data_iterator = haystack(viewFinder)._begin;
            haystack(viewFinder)._end = haystack(viewFinder)._begin;  // Copy might be faster than adding an offset.
            ++entryIt;
            // Iterate over the journal entries.
            for (;entryIt != end(_journalEntries(*it)._journalNodes, Standard()); ++entryIt)
            {
    //            std::cerr << "Processing entry: " << entryIt - begin(_journalEntries(*it)._journalNodes, Standard()) << std::endl;

                switch(entryIt->segmentSource)
                {  // Check for a deletion.
                    case SOURCE_ORIGINAL:
                        haystack(viewFinder)._end += (overlapSize << 1);
    //                    std::cout << "Substring Original: " << infix(*it, haystack(viewFinder)._begin, haystack(viewFinder)._end) << std::endl;
                        break;
                    case SOURCE_PATCH:
                        haystack(viewFinder)._end += entryIt->length + (overlapSize << 1);
    //                    std::cout << "Substring Patch: " << infix(*it, haystack(viewFinder)._begin, haystack(viewFinder)._end) << std::endl;
                        break;
                    default:
                        SEQAN_ASSERT_FAIL("Unknown segment source!");
                }
                while(find(viewFinder, pattern))
                {
                    SEQAN_OMP_PRAGMA(critical)
                    findContext(it - begin(journalSet, Standard()), position(viewFinder) + position(haystack(viewFinder)._begin));
                }

                haystack(viewFinder)._begin += entryIt->length;
                hostIterator(viewFinder).data_iterator = haystack(viewFinder)._begin;  // Need to update the data_iterator manually.
                haystack(viewFinder)._end = haystack(viewFinder)._begin;
            }
            clear(viewFinder);
         }
    }
    std::cout << "Time: " << sysTime() - timeStart << " s.\n";

}

template <typename TSetSpec, typename TNeedle, typename TSpec>
inline void _doFindPattern(FindContext_ & findContext,
                           StringSet<Dna5String, Owner<TSetSpec> > & sequenceSet,
                           Pattern<TNeedle, TSpec> & pattern,
                           FindOptions const & findOptions)
{
    typedef StringSet<Dna5String, Owner<TSetSpec> > TStringSet;
    typedef typename Value<TStringSet>::Type TStringSetValue;
    typedef typename Iterator<TStringSet, Standard>::Type TStringSetIterator;

    // Start with the serial version here.
    double timeStart = sysTime();

    omp_set_num_threads(findOptions.numThreads);

    Splitter<TStringSetIterator> setSplitter(begin(sequenceSet, Standard()), end(sequenceSet, Standard()), Parallel());
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned job = 0; job < length(setSplitter); ++job)
    {
        printf("Thread %ld of %ld threads.\n", job, length(setSplitter));
        // Variant one: We parse over each sequence of the journaled string and search for the pattern.
        for (TStringSetIterator it = setSplitter[job]; it != setSplitter[job+1]; ++it)
        {
            for (unsigned k = 0; k < 50;++k)
            {
                Finder<TStringSetValue> finder;
                setContainer(finder, *it);
                while(find(finder, pattern))
                {
                    SEQAN_OMP_PRAGMA(critical)
                    findContext(it - begin(sequenceSet, Standard()), position(finder));
                }
                clear(finder);
            }
        }
    }
    std::cout << "Time: " << sysTime() - timeStart << " s.\n";
}

//template <typename TSequenceSet, typename TPattern>
//inline int _findPattern(FindContext_ & findContext,
//                        TSequenceSet & seqSet,
//                        TPattern  & patternSeq,
//                        FindOptions const & findOptions)
//{
//    if (findOptions.method == "simple")
//    {
//        Pattern<Dna5String, Simple> searchPattern(patternSeq);
//        _doFindPattern(findContext, seqSet, searchPattern, findOptions);
//    }
//    else if (findOptions.method == "horspool")
//    {
//        Pattern<Dna5String, Horspool> searchPattern(patternSeq);
//        _doFindPattern(findContext, seqSet, searchPattern, findOptions);
//    }
//    else if (findOptions.method == "shift-and")
//    {
//        Pattern<Dna5String, ShiftAnd> searchPattern(patternSeq);
////        _doFindPattern(findContext, seqSet, searchPattern, findOptions);
//    }
//    else if (findOptions.method == "shift-or")
//    {
//        Pattern<Dna5String, ShiftOr> searchPattern(patternSeq);
////        _doFindPattern(findContext, seqSet, searchPattern, findOptions);
//    }
//    return 0;
//}


template <typename TDeltaMap, typename TStringTreeSpec, typename TNeedle, typename TSearchSpec>
int _findPatternOnline(JournaledStringTree<TDeltaMap, TStringTreeSpec> & stringTree,
                       TNeedle const & needle,
                       FindOptions const & findOptions,
                       ResultWriter_ & writer,
                       TSearchSpec const & /*searchSpec*/)
{
    typedef JournaledStringTree<TDeltaMap, TStringTreeSpec> TStringTree;
    typedef Pattern<TNeedle, TSearchSpec> TPattern;
    typedef Finder_<TStringTree, TPattern, Jst<> > TFinder;

    // Now we can process only blocks of variants.
    double timeBlockAll = sysTime();
    HitCollector hitCollector(length(stringSet(stringTree)));

    reinit(stringTree);

    TPattern pattern;
    setHost(pattern, needle);
    TFinder finder(stringTree);

    if (findOptions.numThreads > 1)
        find(finder, pattern, hitCollector, findOptions.k, Parallel());
    else
        find(finder, pattern, hitCollector, findOptions.k, Serial());
    std::cout << "Time for find: " << sysTime() - timeBlockAll << " s." << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Dump matches.

    CharString patternId = "Pattern";
    if (!writer.writeFormatted(patternId, needle, hitCollector.hitPositionSet))
    {
        std::cerr << "Cannot read output file <" << findOptions.outputFile<< ">!" << std::endl;
        return JSeqTools::FILE_READ_ERROR; // TODO(rmaerker): Change to throw Exception.
    }

    return 0;

//    for (unsigned needleId = 0; needleId != length(needleStore); ++needleId)
//    {
//        resize(collectorStore[needleId]._vpBlock, length(collectorStore[needleId].hitPositionSet), 0, Exact());
//
//        TPattern pattern;
//        appendValue(patternStore, pattern);
//        setHost(patternStore[needleId], needleStore[needleId]);
//
//        TTraverser traverser(stringTree, length(needleStore[needleId]) - findOptions.k);
//        appendValue(traverserStore, traverser);
//
//        TFinder finder2(stringTree);
//        appendValue(finderStore, finder2);
//        _init(finderStore[needleId], patternStore[needleId], findOptions.k, typename ErrorsSupported<TFinderFunctor>::Type());
//        finderStore[needleId]._needReinit = false;
//    }
//
//    unsigned chunkSize = length(keys(container(stringTree)));
//    if (findOptions.chunkSize != -1)
//        chunkSize = findOptions.chunkSize;
//
//    int refLength = length(host(journalData(stringTree)));
//    // Load in blocks and then process patterns per block.
//    unsigned varEnd = length(keys(container(stringTree)));
//    unsigned lastBlock;
//    if (chunkSize > 0)
//        lastBlock = std::ceil(static_cast<double>(varEnd) / static_cast<double>(chunkSize));
//    else
//        lastBlock = 1;

//    for (unsigned blockId = 0; blockId < lastBlock; ++blockId)
//    {
//        double timeBlockLoad = sysTime();
//        unsigned blockBegin = blockId * chunkSize;
//        unsigned blockEnd = _min(length(keys(container(stringTree))), (blockId + 1) * chunkSize);
//
//        // We built the journal for the current block here.
//        setBlockBegin(stringTree, blockBegin);
//        setBlockEnd(stringTree, blockEnd);
//        reinitJournal(stringTree, Parallel());
//
//        std::cerr << "Time for block generating: " << sysTime() - timeBlockLoad << " s" << std::endl;
//
//        double timeBlockSearch = sysTime();
//        // We now initialize the block we are currently running.
//        TPosition hostEndPos = length(host(journalData(stringTree)));
//        if (blockId + 1 != lastBlock)
//            hostEndPos = keys(container(stringTree))[blockEnd];
//
//        // Go over each pattern and search within the current block.
//        for (unsigned j = 0; j < length(finderStore); ++j)
//        {
//            if (blockId == 0) // First initialization.
//            {
//                _initSegment(traverserStore[j], begin(container(stringTree), Rooted()) + blockBegin,
//                             begin(container(stringTree), Rooted()) + blockEnd, 0, hostEndPos);
//            }
//            else  // Reinit the master and the branch node end.
//            {
//                traverserStore[j]._masterItEnd = begin(host(journalData(stringTree)), Rooted()) + hostEndPos;
//                traverserStore[j]._branchNodeItEnd = begin(container(stringTree), Rooted()) + blockEnd;
//            }
//
//            traverse(traverserStore[j], finderStore[j], collectorStore[j]);
//            // TODO(rmaerker): Parallelize!
//            // We load store the current virtual offset for each sequence after we processed the current block.
//            // Note we clear the journal string so that only variants of the current block need to be read.
//            for (unsigned i = 0; i < length(collectorStore[j]._vpBlock); ++i)
//                collectorStore[j]._vpBlock[i] += static_cast<int>(length(journalData(stringTree)[i])) - refLength;
//        }
//        std::cerr << "Time for block searching: " << sysTime() - timeBlockSearch << " s" << std::endl;
//    }
//
//    std::cerr << "Time for all blocks: " << sysTime() - timeBlockAll << " s" << std::endl;

//    timeBlock = sysTime();
//    setBlockBegin(stringTree, 0);
//    setBlockEnd(stringTree, varEnd);
//    reinitJournal(stringTree, Parallel());
//
//    for (unsigned i = 0; i < length(needleStore); ++i)
//    {
//
//        HitCollector hitCol(length(stringSet(stringTree));
//        TPattern pattern;
//        setHost(pattern, needleStore[i]);
//        TFinder finder2(stringTree);
//        find(finder2, pattern, hitCol);
//
//    for (unsigned m = 0; m < length(hitCol.hitPositionSet); ++m)
//    {
//        if (!empty(hitCol.hitPositionSet[m]))
//            std::cout << "Hit in seq " << m << ": " <<  hitCol.hitPositionSet[m][0] << std::endl;
//    }
//    }
//    std::cerr << "Time for all: " << sysTime() - timeBlock << " s" << std::endl;

//    for (unsigned m = 0; m < length(collectorStore[0]._vpBlock); ++m)
//        std::cout << "Seq Length: " << m << ": " << refLength + collectorStore[0]._vpBlock[m] << " vs. " << length(journalData(stringTree)[m]) << std::endl;

//    double time = sysTime();
//    find(finder2, onlinePattern, hitCollector);
//    std::cerr << "Time: " << sysTime() - time << " s" <<std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Dump matches.

//    CharString patternId = "Pattern";
//    ResultWriter_ writer(findOptions.outputFile);
//    for (unsigned i = 0; i< length(collectorStore); ++i)
//    {
//        if (!writer.writeFormatted(patternId, needleStore[i], collectorStore[i].hitPositionSet))
//        {
//            std::cerr << "Cannot read output file <" << findOptions.outputFile<< ">!" << std::endl;
//            return JSeqTools::FILE_READ_ERROR; // TODO(rmaerker): Change to throw Exception.
//        }
//    }
//    return 0;
}

// ----------------------------------------------------------------------------
// Function _findPatternCompare
// ----------------------------------------------------------------------------

template <typename TReadStore, typename TMethod>
inline void
_findPatternCompare(TReadStore const & readStore, FindOptions const & findOptions, TMethod const & /*tag*/)
{
    typedef typename Value<TReadStore>::Type TRead;
    typedef typename Value<TRead>::Type TValue;
    typedef String<TValue, Alloc<> > THostSeq;


    typedef Pattern<THostSeq, TMethod> TPattern;
    typedef Finder<THostSeq> TFinder;

    double timeRefSearchAll = sysTime();
    // Now we need to load the ref seq into blocks.
    std::ifstream refStream;
    refStream.open(toCString(findOptions.seqFastaFile), std::ios_base::in);
    if (!refStream.good())
    {
        std::cerr << "Can't read " << findOptions.seqFastaFile << "!" << std::endl;
        return;
    }

    RecordReader<std::ifstream, SinglePass<> > reader(refStream);
    StringSet<THostSeq> refIdStore;
    StringSet<THostSeq> refStore;
    // Resize the sets.
    resize(refStore, findOptions.blockSize, Exact());
    resize(refIdStore, findOptions.blockSize, Exact());
    ResultWriter_ writer(findOptions.outputFile);
    StringSet<String<unsigned> > hits;

    unsigned counter = 0;
    while (!atEnd(reader))
    {
        double timeLoadBlock = sysTime();
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // First load the current block with sequences.
        while(!atEnd(reader) && counter != findOptions.blockSize)
        {
            if (readRecord(refIdStore[counter], refStore[counter], reader, Fasta()) != 0)
            {
                std::cerr << "Error reading: " << refIdStore[counter] << " at " << counter << std::endl;
                return;
            }
            ++counter;
        }
        std::cout << "Time to load block: " << sysTime() - timeLoadBlock << "s."<< std::endl;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Second run sequentially over each sequence and find all pattern.

        double timeProcessBlock = sysTime();
        // Check each pattern.
        for (unsigned j = 0; j < length(readStore); ++j)
        {
            TPattern pattern(readStore[j]);
            resize(hits, findOptions.blockSize, Exact());
            for (unsigned i = 0; i < counter; ++i)
            {
                double timeSingleRefSearch = sysTime();
                TFinder finder(refStore[i]);
                _findPatternForComparison(hits[i], finder, pattern, length(readStore[j]) - findOptions.k, findOptions,
                                          TMethod());
                std::cout << "Time single ref search: " << sysTime() - timeSingleRefSearch << "s." <<std::endl;
            }
            writer.writeFormatted("Pattern", readStore[j], hits);
            clear(hits);
        }
        counter = 0;
        std::cout << "Time to process block: " << sysTime() - timeProcessBlock << "s."<< std::endl;
    }

    std::cout << "Time all ref search: " << sysTime() - timeRefSearchAll << "s." <<std::endl;
}

// ----------------------------------------------------------------------------
// Function _findPattern
// ----------------------------------------------------------------------------

template <typename TAlphabetSeq, typename TMethod>
int _findPattern(FindOptions const & findOptions,
                 TAlphabetSeq const & /*alphabet*/,
                 TMethod const & /*method*/)
{
//    typedef DeltaStore<size_t, TAlphabetSeq> TDeltaStore;
    typedef DeltaMap<unsigned, TAlphabetSeq> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TStringTree;
    typedef typename Host<TStringTree>::Type THost;

    typedef String<TAlphabetSeq, Alloc<> > TNeedle;
//    typedef Pattern<TNeedle, TMethod> TPattern;
//    typedef typename GetStringSet<TStringTree>::Type TJournalSet;  // Only used for the checking.

    double totalTime = sysTime();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Read the pattern file.

    typedef StringSet<TNeedle, Owner<> > TReadStore;
    typedef StringSet<CharString, Owner<ConcatDirect<> > > TReadIdStore;

    std::ifstream patternFile;
    patternFile.open(toCString(findOptions.filePattern), std::ios_base::in);
    if (!patternFile.good())
    {
        std::cerr << "Cannot read input file <" << findOptions.filePattern << ">!" << std::endl;
        return JSeqTools::FILE_READ_ERROR;
    }

    TReadIdStore readIdStore;
    TReadStore readStore;
    RecordReader<std::ifstream, SinglePass<> > patternReader(patternFile);

    while(!atEnd(patternReader))
    {
        TNeedle tmpRead;
        CharString tmpId;
        int res = readRecord(tmpId, tmpRead, patternReader, Fasta());
        if (res)
        {
            std::cout << "Cannot read pattern" << std::endl;
            return -1;
        }
        appendValue(readIdStore, tmpId, Generous());
        appendValue(readStore, tmpRead, Generous());
    }
    patternFile.close();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reference Search.

    if (findOptions.compare)
    {
        _findPatternCompare(readStore, findOptions, TMethod());
        return 0;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Read the reference file.

    THost ref;
    CharString refId;
    _loadSequenceFasta(refId, ref, findOptions.referenceFile);


    // Differ two things:
    // We need to change the behavior of loading the delta data.

    // TODO(rmaerker): Read from vcf file directly!

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Read delta file.

    double timeReadDelta = sysTime();

    String<CharString> nameStore;
    GdfHeader<> gdfHeader;
    gdfHeader._nameStorePtr = &nameStore;

    TDeltaMap deltaMap;
    readJSeqFile(deltaMap, gdfHeader, refId, findOptions.jseqFile);

    TStringTree stringTree(ref, deltaMap);
    if (findOptions.chunkSize != -1)
        setBlockSize(stringTree, findOptions.chunkSize);
    std::cout << "Time for reading delta map: " << sysTime() - timeReadDelta << " s." << std::endl;

    // Set the threads for parallel processing.
    if (findOptions.numThreads > 1)
        omp_set_num_threads(findOptions.numThreads);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // If Checking is enabled we switch to checking mode.

    if (findOptions.checking)
    {
        _checkCorrectAlgorithm(stringTree, findOptions, TMethod());
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Select the search method.
    ResultWriter_ writer(findOptions.outputFile);
    for (unsigned readId = 0; readId < length(readStore); ++readId)
    {
        std::cout << "## Search read " << readId << ":" << std::endl;
        _findPatternOnline(stringTree, readStore[readId], findOptions, writer, TMethod());
    }
    std::cout << "Total time: " << sysTime() - totalTime << "s." << std::endl;
    return 0;

//    RecordReader<std::ifstream, SinglePass<> > reader(patternFile);
//
//    if (findOptions.checking)
//    {
//        typedef typename Value<TJournalSet>::Type TJournalString;
//        typedef typename Iterator<TJournalString>::Type TJournalIter;
//        typedef Finder<TJournalString> TFinderCheck;
//
//        TJournalSet journalSet;
//        setHost(journalSet, ref);
//        adaptTo(journalSet, deltaMap);
//
//        while (!atEnd(reader))
//        {
//            TNeedle needle;
//            CharString seqId;
//            if (readRecord(seqId, needle, reader, Fasta()) != 0)
//                return JSeqTools::FILE_READ_ERROR;
//
//            double timeStart = sysTime();
//            outputStream << "Pattern: " << seqId << '\n';
//            for (unsigned i = 0; i < length(journalSet); ++i)
//            {
//                TFinderCheck finderCheck(journalSet[i]);
//                TPattern pattern(needle);
//                outputStream << jseqHeader._nameStore[i];
//                while(find(finderCheck, pattern))
//                {
//                    outputStream << '\t' << position(finderCheck);
//                }
//                outputStream << '\n';
//            }
//            outputStream << '\n';
//            std::cerr << "Time: "<< sysTime() - timeStart << " s";
//        }
//        return 0;
//    }


    return 0;

}

//template <typename TAlphabet, typename TMethod>
//int _findPatternJournaled(FindOptions const & findOptions,
//                          TAlphabet const & /*alphabet*/,
//                          TMethod const & /*method*/,
//                          StringTreeSparse const & /*tag*/)
//{
//    typedef String<TAlphabet, Alloc<> > THost;
//    typedef JournaledStringTree<THost, StringTreeSparse> TStringTree;
//    typedef typename VariantData<TStringTree>::Type TDeltaMap;
//
////    typedef String<TAlphabet, Journaled<Alloc<>, SortedArray > > TJournalString;
////    typedef typename Host<TJournalString>::Type THostSequence;
////    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournalSet;
//
//    typedef String<TAlphabet, Alloc<> > TNeedle;
//    typedef Pattern<TNeedle, TMethod> TPattern;
//
//    THost ref;
//    CharString refId;
//    _loadSequenceFasta(refId, ref, findOptions.referenceFile);
//
//    JSeqHeader jseqHeader;
//
//    TDeltaMap deltaMap;
//
//    readJSeqFile(deltaMap, jseqHeader, refId, findOptions.jseqFile);
//
//    // Read the patten file!
//    std::ifstream patternFile;
//    patternFile.open(toCString(findOptions.filePattern), std::ios_base::in);
//    if (!patternFile.good())
//    {
//        std::cerr << "Cannot read input file <" << findOptions.filePattern << ">!" << std::endl;
//        return JSeqTools::FILE_READ_ERROR;
//    }
//
//    std::ofstream outputStream;
//    outputStream.open(toCString(findOptions.outputFile), std::ios_base::out);
//    if (!outputStream.good())
//    {
//        std::cerr << "Cannot read input file <" << findOptions.outputFile<< ">!" << std::endl;
//        return JSeqTools::FILE_READ_ERROR; // TODO(rmaerker): Change to throw Exception.
//    }
//
//    RecordReader<std::ifstream, SinglePass<> > reader(patternFile);
//    HitCollector hitCollector(jseqHeader._nameStore);
//    TStringTree stringTree(ref, deltaMap);
//
//    // Search all patterns from left to right.
//    while (!atEnd(reader))
//    {
//        hitCollector.reinit();
//
//        TNeedle needle;
//        CharString seqId;
//        if (readRecord(seqId, needle, reader, Fasta()) != 0)
//            return JSeqTools::FILE_READ_ERROR;
//
//        TPattern pattern;
//        setHost(pattern, needle);
//        _patternInit(pattern);
//
//        Finder_<TStringTree, TPattern, Jst<> > finder2(stringTree, pattern);
//
//        double time = sysTime();
//        find(finder2, pattern, hitCollector);
//        std::cerr << "Time: " << sysTime() - time << " s" <<std::endl;
//        std::cout << "# Search pattern: " << seqId << " = " << needle << "\n" << std::flush;
//        hitCollector.writeHits(outputStream, seqId);
//        std::cout << "\n" << std::flush;
//    }
//    patternFile.close();
//    outputStream.close();
//
////    typedef typename Value<TJournalSet>::Type TJournalString;
////    typedef typename Iterator<TJournalString>::Type TJournalIter;
////    typedef typename Iterator<THost>::Type THostIt;
////
////    TJournalIter it = begin(journalSet[100]);
////    TJournalIter itEnd = end(journalSet[100]);
////
////    double itTime = sysTime();
////    TAlphabet alph;
////    for (; it != itEnd; ++it)
////        alph = getValue(it);
////    std::cerr << "Alph: " << alph << std::endl;
////    std::cerr << "Time It: " << sysTime() - itTime << std::endl;
////
////    THostIt itHost = begin(host(journalSet));
////    THostIt itHostEnd = end(host(journalSet));
////
////    itTime = sysTime();
////    for (; itHost != itHostEnd; ++itHost)
////        alph = getValue(itHost);
////    std::cerr << "Alph: " << alph << std::endl;
////    std::cerr << "Time It: " << sysTime() - itTime << std::endl;
////
////
////    std::cerr << "Find in multiple genomes not parallel!" << std::endl;
////    FindContext_ findContext;
////    bool res = _findPattern(findContext, journalSet, needle, findOptions);
////    findContext.print();
//
//
//    return 0;
//}

// ----------------------------------------------------------------------------
// Function _findPattern
// ----------------------------------------------------------------------------

//int _findPatternVcf(FindOptions const & findOptions)
//{
//    typedef String<Dna5, Journaled<Alloc<>, SortedArray > > TJournalString;
//    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournalSet;
//
//    Dna5String reference;
//    CharString refId;
//
//    int res = readSequence(refId, reference, findOptions.referenceFile);
//    if (res)
//        return res;
//
//    Dna5String patternSeq;
//    CharString patternId;
//    res = readSequence(patternId, patternSeq, findOptions.filePattern);
//    if (res != 0)
//    {
//        std::cerr << "Can't read file <" << findOptions.filePattern << "!"<< std::endl;
//        return 1;
//    }
//
//    VcfHeader vcfHeader;
//    DeltaHistoryTree<Dna5String> deltaHistoryTree(length(patternSeq), 1092, reference);
//    double timeReadingVcf = sysTime();
//    _convertVcfToDeltaHistoryTree(deltaHistoryTree, vcfHeader, findOptions);
//    std::cout << "Time for parsing vcf: " << sysTime() - timeReadingVcf << " s" << std::endl;
//
//
////
////    patternSeq = infix(globalReference(container(deltaHistoryTree)), 34123456,34123476);
////
////    double timeStart = sysTime();
////    find(deltaHistoryTree, patternSeq);
////    std::cout << "Time for search: " << sysTime() - timeStart << " s" << std::endl;
//////    insert(journalSet[0], 12345, patternSeq);
////
////    // Have JournalSet and the reference. Now we can search for the pattern.
////    FindContext_ findContext;
////    res = _findPattern(findContext, journalSet, patternSeq, findOptions);
////    if (res)
////        return res;
////
////    for (unsigned i =0; i < length(findContext.hitString); ++i)
////    {
////        std::cout << "Sequence: " << findContext.hitString[i].i1 << "\t";
////        std::cout << "Position :" << findContext.hitString[i].i2 << std::endl;
////    }
//
//
//    return 0;
//}

int _findPatternFasta(FindOptions const & /*findOptions*/)
{

//    Dna5String reference;
//    CharString refId;
//
//    int res = readSequence(refId, reference, findOptions.referenceFile);
//    if (res)
//        return res;
//
//    std::ifstream seqStream(toCString(findOptions.cglFile), std::ios_base::in);
//    if (!seqStream.good())
//    {
//        std::cerr << "ERROR: Could not open the file.\n";
//        return 1;
//    }
//
//    Dna5String patternSeq = infix(reference, 34123456,34123466);
//    RecordReader<std::ifstream, SinglePass<> > recReader(seqStream);
//
//    TStringSet sequenceSet;
//    TIdStringSet idSet;
//
//    if (findOptions.chunkSize == -1)
//    {
//        while(!atEnd(recReader))
//        {
//            Dna5String tmpSeq;
//            CharString tmpId;
//            if (readRecord(tmpId, tmpSeq, recReader, Fasta()) != 0)
//            {
//                std::cerr << "Error when reading the fasta file." << std::endl;
//                return 1;
//            }
//            appendValue(sequenceSet, tmpSeq, Generous());
//            appendValue(idSet, tmpId, Generous());
//        }
//
//        FindContext_ findContext;
//        _findPattern(findContext, sequenceSet, patternSeq, findOptions);
//
//        for (unsigned i =0; i < length(findContext.hitString); ++i)
//        {
//            std::cout << "Sequence: " << findContext.hitString[i].i1 << "\t";
//            std::cout << "Position :" << findContext.hitString[i].i2 << std::endl;
//        }
//        return 0;
//    }
//
//    resize(sequenceSet, findOptions.chunkSize, Exact());
//    resize(idSet, findOptions.chunkSize, Exact());
//
//    unsigned currSeqBlock = 0;
//    while (!atEnd(recReader))
//    {
//        unsigned counter = 1;
//        // We read the sequences one after another.
//        for (unsigned i = 0; i < length(sequenceSet); ++i, ++counter)
//        {
//            Dna5String tmpSeq;
//            CharString tmpId;
//            if (readRecord(tmpId, tmpSeq, recReader, Fasta()) != 0)
//            {
//                std::cerr << "Error when reading the fasta file." << std::endl;
//                return 1;
//            }
//            sequenceSet[i] = tmpSeq;
//            idSet[i] = tmpId;
//
//            if (atEnd(recReader))
//            {
//                resize(sequenceSet, counter, Exact());
//                resize(idSet, counter, Exact());
//                break;
//            }
//        }
//        FindContext_ findContext;
//        _findPattern(findContext, sequenceSet, patternSeq, findOptions);
//
//        for (unsigned i =0; i < length(findContext.hitString); ++i)
//        {
//            std::cout << "Sequence: " << findContext.hitString[i].i1 + currSeqBlock << "\t";
//            std::cout << "Position :" << findContext.hitString[i].i2 << std::endl;
//        }
//        currSeqBlock+=length(sequenceSet);
//    }

    return 0;
}

}

#endif // EXTRAS_TOOLS_JSEQ_TOOLS_JSEQ_TOOLS_FIND_H_
