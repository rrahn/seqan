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
// This tool implements online-pattern search over multiple genomes.
// ==========================================================================

//#ifdef SEQAN_ENABLE_DEBUG
//#undef SEQAN_ENABLE_DEBUG
//#endif
//#define SEQAN_ENABLE_DEBUG 0
//#define PROFILE_DATA_PARALLEL_INTERN
//#define DEBUG_DATA_PARALLEL_2
//#define DEBUG_DATA_PARALLEL

#define PROFILE_JST_INTERN

#include "jst_bench_facade_header.h"
#include "jst_bench_find.h"   // Finding pattern in set.

using namespace seqan;

// ----------------------------------------------------------------------------
// Function readCheckOptions
// ----------------------------------------------------------------------------

template <typename TParser>
inline ArgumentParser::ParseResult
readCheckOptions(FindCheckerOptions & options, TParser const & parser)
{
    if (isSet(parser, "selftest"))
        options.checking = true;
    if (isSet(parser, "window-length"))
    {
        getOptionValue(options.wlRange.i1, parser, "window-length", 0);
        getOptionValue(options.wlRange.i2, parser, "window-length", 1);
    }
//    if (isSet(parser, "errors"))
//        getOptionValue(options.errors, parser, "error-rate");
    if (isSet(parser, "rng"))
        getOptionValue(options.seed, parser, "rng");
    if (isSet(parser, "cycles"))
        getOptionValue(options.patternCount, parser, "cycles");
    return ArgumentParser::PARSE_OK;
}

template <typename TParser>
inline ArgumentParser::ParseResult
readCompareOptions(FindComparisonOptions & options, TParser const & parser)
{
    if (isSet(parser, "sequential"))
        options.compare = true;


    if (!isSet(parser, "si"))
    {
        std::cerr << "Invalid multi fasta file" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    getOptionValue(options.seqFastaFile, parser, "si");
    getOptionValue(options.blockSize, parser, "b");
    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function readFinderOptions()
// ----------------------------------------------------------------------------

inline ArgumentParser::ParseResult
readFinderOptions(FindOptions & options, ArgumentParser const & parser)
{

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Read general program options.

    // Read the main options of the program.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.jseqFile, parser, "input");
    getOptionValue(options.referenceFile, parser, "reference");
    getOptionValue(options.filePattern, parser, "pattern");
    getOptionValue(options.outputFile, parser, "output");
    getOptionValue(options.numThreads, parser, "threads");
    // Parallelize with OpenMp
    omp_set_num_threads(options.numThreads);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Online search program options.

    getOptionValue(options.method, parser, "algorithm");
    CharString alphabet;
    getOptionValue(alphabet, parser, "pattern-alphabet");
    if (alphabet == "dna")
        options.alphabetId = 0;
    else if(alphabet == "dna5")
        options.alphabetId = 1;

    if (isSet(parser, "errors"))
    {
        int tmp;
        getOptionValue(tmp, parser, "errors");
        options.k = -tmp;
    }


    if (isSet(parser, "cs"))
        getOptionValue(options.chunkSize, parser, "cs");

    if (isSet(parser, "selftest"))
        readCheckOptions(options, parser);

    if (isSet(parser, "sequential"))
        readCompareOptions(options, parser);

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function _setupAndParseArguments()
// ----------------------------------------------------------------------------

template <typename TArgCount, typename TArgValues>
ArgumentParser::ParseResult _setupAndParseArguments(FindOptions & options,
                             TArgCount argc,
                             TArgValues const & argv)
{
    ArgumentParser parser("jst_bench");
    // Set short description, version, and date.
    setShortDescription(parser, "Journaled String Tree Benchmark Tool.");
    setVersion(parser, "1.1");
    setDate(parser, "April 2014");

    // Define usage line and long description.
    addUsageLine(parser," \"[\\fIOPTIONS\\fP] \\fB-i\\fP <IN> \\fB-r\\fP <REF> \\fB-p\\fP <PATTERN> \\fB-o\\fP <OUT>\" ");
    addDescription(parser, "Benchmarking tool to evaluate simultaneous traversal over multiple sequences using the Journaled String Tree.");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section General.

    addSection(parser, "General Options");
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("i", "input", "GDF file.", ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, "input", "gdf");
    setRequired(parser, "input", true);

    addOption(parser, seqan::ArgParseOption("r", "reference", "File with reference sequence.", ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, "reference", "fa fasta");
    setRequired(parser, "reference", true);

    addOption(parser, ArgParseOption("p", "pattern", "Fasta file containing the pattern.", ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, "pattern", "fa fasta");
    setRequired(parser, "pattern", true);

    addOption(parser,ArgParseOption("o", "output", "Output file", ArgParseArgument::STRING));
    setRequired(parser, "output", true);

    addOption(parser, ArgParseOption("t", "threads", "Number of threads used for parallel processing.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "threads", "1");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section Online Search.

    addSection(parser, "JST Search");

    addOption(parser, ArgParseOption("pa", "pattern-alphabet", "Alphabet used for the pattern.", ArgParseArgument::STRING));
    setValidValues(parser, "pattern-alphabet", "dna dna5");
    setDefaultValue(parser, "pattern-alphabet", "dna");

    addOption(parser, ArgParseOption("a", "algorithm", "Algorithm used for searching the pattern.", ArgParseArgument::STRING));
    setValidValues(parser, "algorithm", "simple horspool shift-and shift-or myers");
    setDefaultValue(parser, "algorithm", "horspool");

    addOption(parser, ArgParseOption("e", "errors", "Number of errors allowed for approximate search.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "e", "0");
    setMinValue(parser, "e", "0");
    setMaxValue(parser, "e", "10");

    addOption(parser, ArgParseOption("cs", "chunk-size", "The number of variants processed per cycle.", ArgParseArgument::INTEGER));
    setMinValue(parser, "chunk-size", "10000");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section Reference Search

    addSection(parser, "Sequential Search");
    addOption(parser, ArgParseOption("","sequential", "Flag to indicate sequential search."));
    addOption(parser, ArgParseOption("si", "sequential-input", "Multi-fasta file containing the sequences to process.", ArgParseArgument::INPUTFILE));
    setValidValues(parser, "si", "fa fasta");
    addOption(parser, ArgParseOption("b", "blocks", "Number of sequences processed per cycle.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "b", "100");

//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    // Section Filtration.
//
//    addSection(parser, "Multiple Pattern Search");
//
//    addOption(parser, ArgParseOption("", "filtration", "Enables approximate pattern search with filtration."));
//    setHelpText(parser, "filtration", "Use this option if you want to search many patterns.");
//
//    addOption(parser, ArgParseOption("s", "seed", "Size of the seeds to build the index over the patterns.", ArgParseArgument::INTEGER));
//    // TODO(rmaerker): Add option for qgram interleave
//    // TODO(rmaerker): Add option for Hamming distance or indel search.
//    addOption(parser, ArgParseOption("er", "error-rate", "The error rate to allow for approximate hits.", ArgParseArgument::DOUBLE));
//    setMinValue(parser, "error-rate", "0.0");  // Allow exact search too.
//    setMaxValue(parser, "error-rate", "0.1");  // Maximal error rate is 10%.
//
//    setDefaultValue(parser, "error-rate", "0.02");  // Default error rate is 5 %.
//    addOption(parser, ArgParseOption("", "hamming", "Uses hamming distance for verification, edit distance otherwise."));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section Selftest.

    // TODO(rmaerker): Hide if done with the tool
    addSection(parser, "Selftest");

    addOption(parser, ArgParseOption("", "selftest", "Forces execution of a selftest."));

    addOption(parser, ArgParseOption("wl", "window-length", "The range of the window size.", ArgParseArgument::INTEGER, "BEGIN END", false, 2));
    setMinValue(parser, "window-length", "1");
    setMaxValue(parser, "window-length", "500");

    addOption(parser, ArgParseOption("c", "cycles", "The number of independent runs.", ArgParseArgument::INTEGER));

    addOption(parser, ArgParseOption("","rng", "The seed for the random number generator.", ArgParseArgument::INTEGER));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Example.

    // Add Examples Section.
//    addTextSection(parser, "Examples");
//    addListItem(parser, "\\fBjst_bench\\fP \\fB-i\\fP DELTA.gdf \\fB-r\\fP REF.fa \\fB-p\\fP PATTERN.fa \\fB-o\\fP OUT.txt \\fP\\fB-a \\fP\\fImyers \\fP\\fB-e \\fP\\fI3\\fP",
//                        "Call with \\fITEXT\\fP set to \"text\" with verbose output.");



    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    return readFinderOptions(options, parser);
}

template <typename TAlphabet>
inline int _findPattern(FindOptions const & options,
                        TAlphabet const & alphabet)
{
    if (options.approximate)
    {
//        return _findApproximative
    }


    if (options.method == "simple")
        return _findPattern(options, alphabet, Simple());
    else if (options.method == "horspool")
        return _findPattern(options, alphabet, Horspool());
    else if (options.method == "shift-and")
        return _findPattern(options, alphabet, ShiftAnd());
    else if (options.method == "shift-or")
        return _findPattern(options, alphabet, ShiftOr());
    else if (options.method == "myers")
        return _findPattern(options, alphabet, Myers<FindInfix, True, void>());
    else
        return -1;

    // TODO(rmaerker): Add myers bit vector.
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const** argv)
{

    FindOptions options;
    ArgumentParser::ParseResult res = _setupAndParseArguments(options, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
    {
        std::cerr << "Errors while parsing the input occurred!" << std::endl;
        return res;
    }

// NOTE(rmaerker): This is debug code to extract a reference sequence.
//    typedef StringSet<Dna5String, Owner<> > TReadStore;
//    typedef StringSet<CharString, Owner<ConcatDirect<> > > TReadIdStore;
//
//    std::ifstream file;
//    file.open(toCString(options.referenceFile), std::ios_base::in);
//    if (!file.good())
//    {
//        std::cerr << "Cannot read input file <" << options.referenceFile << ">!" << std::endl;
//        return JSeqTools::FILE_READ_ERROR;
//    }
//
//    TReadIdStore readIdStore;
//    TReadStore readStore;
//    RecordReader<std::ifstream, SinglePass<> > patternReader(file);
//
//    std::ofstream oFile;
//    oFile.open("chr22_GRCh37.fa", std::ios_base::out);
//    if (!oFile.good())
//    {
//        std::cerr << "Cannot write to output file <" << "chr22_GRCh37.fa" << ">!" << std::endl;
//        return JSeqTools::FILE_READ_ERROR;
//    }
//
//    unsigned count = 0;
//    while(!atEnd(patternReader))
//    {
//        Dna5String tmpRead;
//        CharString tmpId;
//        int res = readRecord(tmpId, tmpRead, patternReader, Fasta());
//        if (res)
//        {
//            std::cerr << "Cannot read pattern" << std::endl;
//            return -1;
//        }
//        appendValue(readIdStore, tmpId, Generous());
//        if (count == 21)
//        {
//            writeRecord(oFile, tmpId, tmpRead, Fasta());
//        }
//        ++count;
////        appendValue(readStore, tmpRead, Generous());
//    }
//    file.close();
//    oFile.close();
//    for (unsigned i = 0; i < length(readIdStore); ++i)
//    {
//        std::cout << "Pos: " << i << " " << readIdStore[i] << std::endl;
//    }

    switch (options.alphabetId)
    {
        case 0:
            _findPattern(options, Dna());
            break;
        case 1:
            _findPattern(options, Dna5());
            break;
        default:
            return -1;
    }
}
