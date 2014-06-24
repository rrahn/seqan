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

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

#include "jst_mapper.h"
#include "jst_mapper_map_reads.h"

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class MyFragStoreConfig
// ----------------------------------------------------------------------------

// Copied from razers3.cpp - Check need for this tool.
struct MyFragStoreConfig :
    public FragmentStoreConfig<>
{
    typedef String<Dna5> TContigSeq;
};

// ==========================================================================
// Functions
// ==========================================================================

template <typename TVerfiierSpec>
JstMapperResult
run(JstMapperOptions const & options)
{
    // Later we need a DeltaContigMap<> -> to map between the contig and it's delta file.
    typedef FragmentStore<MyFragStoreConfig>                        TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore                  TReadStore;
    typedef typename TFragmentStore::TContigStore                   TContigStore;
    typedef typename Value<TContigStore>::Type                      TContigStoreElement;
    typedef typename TContigStoreElement::TContigSeq                TContigSeq;
    typedef typename Value<TContigSeq>::Type                        TContigSeqAlphabet;
    typedef typename Position<TContigSeq>::Type                     TContigSeqPosition;

    typedef DeltaMap<TContigSeqPosition, TContigSeqAlphabet>        TDeltaMap;

    // 1) Load contigs, deltas and reads.
    FragmentStore<MyFragStoreConfig> fragStore;
    TDeltaMap deltaMap;

    if (!loadContigs(fragStore, options.inputGenome))
        return LOAD_CONTIGS_FAILED;

    if (!load(deltaMap, options.inputGenomeDelta))
        return LOAD_DELTA_FAILED;

    if (!loadReads(fragStore, options.inputReads1))
        return LOAD_READS_FAILED;

    // 2) Build an index over all reads and a SWIFT pattern over this index.
    mapReads<TVerfiierSpec>(fragStore, deltaMap, options);

    TIndex index(fragStore.readSeqStore);
    TPattern pattern(index);

    // 3) Enumerate all epsilon matches.
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i) {
        TFinder finder(fragStore.contigStore[i].seq);
        while (find(finder, pattern, EPSILON)) {
            // Verify match.
            Finder<TContigSeq> verifyFinder(fragStore.contigStore[i].seq);
            setPosition(verifyFinder, beginPosition(finder));
            Pattern<TReadSeq, HammingSimple> verifyPattern(fragStore.readSeqStore[position(pattern).i1]);
            unsigned readLength = length(fragStore.readSeqStore[position(pattern).i1]);
            int minScore = -static_cast<int>(EPSILON * readLength);
            while (find(verifyFinder, verifyPattern, minScore) && position(verifyFinder) < endPosition(infix(finder))) {
                TAlignedRead match(length(fragStore.alignedReadStore), position(pattern).i1, i,
                                   beginPosition(verifyFinder), endPosition(verifyFinder));
                appendValue(fragStore.alignedReadStore, match);
            }
        }
    }

    // 4) Write out Sam file.
    std::ofstream samFile(argv[3], std::ios_base::out);
    write(samFile, fragStore, Sam())

}

inline JstMapperResult
configureAndRun(JstMapperOptions const & mapperOptions)
{
    switch (mapperOptions.errorModel)
    {
        case VerifierErrorModel::HAMMING_DISTANCE:
            return run<Hamming_>(mapperOptions);
        case VerifierErrorModel::EDIT_DISTANCE:
            return run<EditDistance>(mapperOptions);
        default:
            return UNKNOWN_OPTION_VALUE;
    }
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

void setUpArgumentParser(ArgumentParser & parser)
{
    // Setup ArgumentParser.
    setAppName("jst_mapper");
    // Set short description, version, and date.
    setShortDescription(parser, "Multiple Reference Mapper");
    setVersion(parser, "0.1");
    setDate(parser, "Juni 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIGENOME\\fP \\fIGENOME-DELTA\\fP \\fIREADS\\fP");
    addDescription(parser, "Write me.");

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "GENOME"));
    setValidValues(parser, 0, "fa fasta");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "GENOME-DELTA"));
    setValidValues(parser, 1, "gdf");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "READS"));
    setValidValues(parser, 2, "fasta fa fastq");


    // Main Options.
    addOption(parser, seqan::ArgParseOption("o", "output", "File to write the matches to."));
    setValidValues(parser, "o", "bjst");

    addOption(parser, ArgParseOption("js", "jst-size", "Number of deltas journaled at a time.", ArgParseArgument::INTEGER));
    setMinValue(parser, "js", "10000");

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Filter Options.
    addOption(parser, ArgParseOption("qs", "qGram-size", "Size of the q-grams used for filtering.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "qs", "15");
    setMinValue(parser, "qs", 8);
    setMaxValue(parser, "qs", 30);


    // Verifier Option
    addOption(parser, ArgParseOption("", "hamming", "Enables verification using hamming distance. If not set edit distance is used"));
    addOption(parser, ArgParseOption("er", "error-rate", "Maximal allowed error rate to find matches."), ArgParseArgument::DOUBLE);
    setDefaultValue(parser, "er", "0.0");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBjst_mapper\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");
}

void extractOptions(JstMapperOptions & options, ArgumentParser const & parser)
{
    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getArgumentValue(options.inputGenome, parser, 0);
    getArgumentValue(options.inputGenomeDelta, parser, 1);
    getArgumentValue(options.inputReads1, 2);

    if (isSet(parser, "js"))
        getOptionValue(options.jstSize, parser, "js");
    else
        options.jstSize = 0;

    // Extract Filter Options.
    getOptionValue(options.qGram, parser, "qs");

    // Extract Verifier Options.
    if (isSet(parser, "hamming"))
        options.errorModel = HAMMING_DISTANCE;
    else
        options.errorModel = EDIT_DISTANCE;

    getOptionValue(options.errorRate, parser, "er");
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    JstMapperOptions options;

    setUpArgumentParser(parser);

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;


    configureAndRun(options);

    return 0;
}
