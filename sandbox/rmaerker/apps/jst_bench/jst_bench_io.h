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
// Implements the basic io features for reading/writing data from
// from/to journal sequences.
// ==========================================================================

#ifndef EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_IO_H_
#define EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_IO_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct ReadVariantInformation_;
typedef Tag<ReadVariantInformation_> WithoutLoadingJournalData;

struct ReadVariantInformationAndJournalData_;
typedef Tag<ReadVariantInformationAndJournalData_> WithLoadingJournalData;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readJSeqFile()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TConfig, typename TFileLocation>
inline int readJSeqFile(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
                        GdfHeader<TConfig> & gdfHeader,
                        CharString const & refId,
                        TFileLocation const & filePath)
{
    std::ifstream inputStream;
    inputStream.open(toCString(filePath), std::ios_base::in);
    if (!inputStream.good())
    {
        std::cerr << "Cannot read file <" << filePath << ">!" << std::endl;
        return JSeqTools::FILE_READ_ERROR;
    }

//    RecordReader<std::ifstream, SinglePass<> > reader(inputStream);

    read(deltaMap, gdfHeader, inputStream, Gdf());

    inputStream.close();
    // Simple reference checking:
    if (gdfHeader._refInfos._refId != refId)
        return GdfIO::UNSUPPORTED_REFERNCE_INFORMATION_ERROR;  // TODO(rmaerker): Change the name of this error message.

    return JSeqTools::FILE_READ_OK;
}

// ----------------------------------------------------------------------------
// Function _loadSequenceFasta()
// ----------------------------------------------------------------------------

template <typename TId, typename TSequence>
inline int _loadSequenceFasta(TId & idString,
                              TSequence & sequence,
                              CharString file)
{
    std::ifstream fileStream(toCString(file), std::ios_base::in);
    if (!fileStream.good())
    {
        std::cerr << "Can't read file <" << file << "!"<< std::endl;
        return 1;
    }

    RecordReader<std::ifstream, SinglePass<> > recReader(fileStream);
    int res = readRecord(idString, sequence, recReader, Fasta());
    if (res != 0)
        return JSeqTools::FILE_READ_ERROR;
    return JSeqTools::FILE_READ_OK;
}

}

#endif // EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_IO_H_
