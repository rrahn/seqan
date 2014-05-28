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
// Implements the header for the journal sequence file formats.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_GDF_IO_GDF_IO_HEADER_H_
#define EXTRAS_INCLUDE_SEQAN_GDF_IO_GDF_IO_HEADER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class JSeqHeaderRecord
// ----------------------------------------------------------------------------

struct JSeqHeaderRecord
{
    CharString _key;
    CharString _value;

    JSeqHeaderRecord() : _key(""), _value("")
    {}

    JSeqHeaderRecord(CharString const & key, CharString const & value) : _key(key), _value(value)
    {}
};

struct JSeqHeaderReferenceBlock
{
    unsigned    _refHash;
    CharString  _refId;
    CharString  _refFile;
};

struct JSeqHeaderFileBlock
{
    unsigned    _minorFileId;
    unsigned    _majorFileId;
    unsigned    _blockSize;
    bool        _byteOrder; // true = little endian; false = big endian.
    bool        _snpCompression;  // true if dna or rna, otherwise false.

    JSeqHeaderFileBlock() : _minorFileId(JSeqIO::FILE_VERSION_LITTLE),
                            _majorFileId(JSeqIO::FILE_VERSION_BIG),
                            _blockSize(JSeqFileBlockSize_<void>::VALUE),
                            _byteOrder(SystemsByteOrder::IS_LITTLE_ENDIAN()),
                            _snpCompression(false)
    {}
};


// ----------------------------------------------------------------------------
// Class JSeqHeader
// ----------------------------------------------------------------------------

class JSeqHeader
{
public:
    JSeqHeaderFileBlock _fileInfos;
    JSeqHeaderReferenceBlock _refInfos;

    String<CharString>  _nameStore;  // The names of the individuals
    String<JSeqHeaderRecord> _headerRecords;  // Variable informations.

    JSeqHeader() :
        _fileInfos(),
        _refInfos(),
        _nameStore(),
        _headerRecords()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setKey()
// ----------------------------------------------------------------------------

inline void
setKey(JSeqHeaderRecord & record, CharString const & key)
{
    record._key = key;
}

// ----------------------------------------------------------------------------
// Function getKey()
// ----------------------------------------------------------------------------

inline CharString &
getKey(JSeqHeaderRecord & record)
{
    return record._key;
}

inline CharString const &
getKey(JSeqHeaderRecord const & record)
{
    return record._key;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

inline void
setValue(JSeqHeaderRecord & record, CharString const & value)
{
    record._value = value;
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

inline CharString &
getValue(JSeqHeaderRecord & record)
{
    return record._value;
}

inline CharString const &
getValue(JSeqHeaderRecord const & record)
{
    return record._value;
}

}

#endif // EXTRAS_INCLUDE_SEQAN_GDF_IO_GDF_IO_HEADER_H_
