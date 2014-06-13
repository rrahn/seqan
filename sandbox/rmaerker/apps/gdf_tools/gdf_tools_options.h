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
// Defines the options used for the different operations available for
// the jseq_tools.
// ==========================================================================

#ifndef EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_OPTIONS_H_
#define EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_OPTIONS_H_


namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// The general app options.
struct AppOptions
{
    unsigned verbosity;  // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    // TODO(rmaerker): referenceFile, jseqIndexFile, bjseq file, djseq file
    unsigned numThreads;

    CharString referenceFile;        // We need the reference file (should be an external file!)! It's mandatory!
    CharString jseqFile;
    CharString cglIndexFile;     // This is the input index file which we might want to generate.

    AppOptions() : verbosity(1), numThreads(1)
    {}
};

// --------------------------------------------------------------------------
// Class Converter
// --------------------------------------------------------------------------

struct ConverterOptions
{
    unsigned verbosity;  // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    unsigned numThreads;

    CharString inputFile;
    CharString outputFile;

    unsigned refAlphabet;
    unsigned varAlphabet;

    // VCF Options:
    CharString vcfReferenceFile;
    bool readGenotype;
    bool includeReference;
    String<unsigned> haplotypes;


    unsigned numIndividuals;

    bool suppressSVs;

    ConverterOptions() : verbosity(1), numThreads(1), refAlphabet(), varAlphabet(), readGenotype(false), includeReference(false), numIndividuals(1), suppressSVs(false)
    {}
};

// ----------------------------------------------------------------------------
// Struct CompressorOptions
// ----------------------------------------------------------------------------

struct JSeqCompressOptions
{
    CharString refSeqFile;
    CharString refIndexFile;
    CharString multiFastaFile;
    CharString outputFile;

    unsigned threads;

    unsigned seedLength;
    unsigned maxGapSize;

    JSeqCompressOptions() : threads(), seedLength(30), maxGapSize(10000000)
    {}
};


// --------------------------------------------------------------------------
// Class ViewerOptions
// --------------------------------------------------------------------------

struct ViewOptions : AppOptions
{
    unsigned streamFormat;  // 0 = no formatting, 1 - FASTA
    String<int> seqPositions;      // Print the file with these position identifiers.
    StringSet<CharString> seqNames;   // Print the file with these name identifiers.
    CharString outputFile;      // Write the output to, default: stdout


    ViewOptions()
        : AppOptions(), streamFormat(0), seqPositions(), seqNames()
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}

#endif // EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_OPTIONS_H_