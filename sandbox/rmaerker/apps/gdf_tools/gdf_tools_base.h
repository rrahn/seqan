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
// Facade header for the jseq_tools.
// ==========================================================================

#ifndef EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_BASE_H_
#define EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_BASE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct JSeqTools
{
    enum Command
    {
        COMMAND_DISPLAY,
        COMMAND_FIND,
        COMMAND_COMPRESS,
        COMMAND_UNKNOWN = -1
    };

    enum FileHandling
    {
        FILE_FORMAT_NOT_SUPPORTED_ERROR,
        FILE_READ_ERROR,
        FILE_READ_OK
    };
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// TODO(rmaerker): Dummy function: Need to be changed later.
template <typename TResult, typename TSource>
inline TResult _readInteger(TSource const & /*val*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function _printJournalStat()
// ----------------------------------------------------------------------------

template <typename TJournalString>
inline void _printJournalStat(TJournalString const & journal)
{
    typedef typename JournalType<TJournalString const>::Type TJournalEntries;
    typedef typename Value<TJournalEntries>::Type TJournalNode;

    std::cerr << "Number of journal entries: " << length(journal._journalEntries._journalNodes)  << std::endl;
    std::cerr << "Lenght of Insertionbuffer: " << length(journal._insertionBuffer) << std::endl;
    std::cerr << "sizeof(TJournal()) " << sizeof(TJournalNode) << std::endl;
    unsigned memUsage = sizeof(TJournalNode) * length(journal._journalEntries._journalNodes) + length(journal._insertionBuffer);
    std::cerr << "Memory usage: " <<  memUsage + sizeof(TJournalString) << " Bytes." << std::endl;
    std::cerr << "Memory ref: " <<  length(host(journal)) << " Bytes." << std::endl;
}

// ----------------------------------------------------------------------------
// Function readGeneralOptions()
// ----------------------------------------------------------------------------

inline void
readGeneralOptions(AppOptions & options, ArgumentParser const & parser)
{
    // Read the main options of the program.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    getOptionValue(options.jseqFile, parser, "compressed-library");
    // TODO(rmaerker): Check correct input file format!
    getOptionValue(options.referenceFile, parser, "reference");
}


}

#endif // EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_BASE_H_
