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
// Implements strategy to adapt a vcf file to journaled data.
// ==========================================================================

#ifndef EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_ADAPT_VCF_H_
#define EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_ADAPT_VCF_H_

#include <seqan/vcf_io.h>
#include "gdf_tools_facade_header.h"

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename T>
struct VcfRecordTranslator{};

template <typename TValue, typename TAlphabet>
struct VcfRecordTranslator<DeltaMap<TValue, TAlphabet> >
{

    typedef DeltaMap<TValue, TAlphabet> TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type TSnp;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;

    TDeltaMap * _variantStorePtr;

    VcfRecordTranslator() : _variantStorePtr(0)
    {}

    VcfRecordTranslator(TDeltaMap & variantStore) : _variantStorePtr(&variantStore)
    {}

    template <typename TVariantSet, typename TCoverageSet>
    inline void operator()(TVariantSet const & varSet, TCoverageSet const & coverageSet, bool force = false)
    {
        SEQAN_ASSERT_EQ(length(varSet), length(coverageSet));

        for (unsigned j = 0; j < length(coverageSet); ++j)
        {
            if (!force && testAllZeros(coverageSet[j]))
                continue;

            switch(varSet[j].variantType)
            {
                case DeltaType::DELTA_TYPE_SNP:
                    insert(*_variantStorePtr, varSet[j].refPos, static_cast<TSnp>(varSet[j].seqBuffer[0]), coverageSet[j]);
                    break;
                case DeltaType::DELTA_TYPE_INS:
                    insert(*_variantStorePtr, varSet[j].refPos, varSet[j].seqBuffer, coverageSet[j]);
                    break;
                case DeltaType::DELTA_TYPE_DEL:
                    insert(*_variantStorePtr, varSet[j].refPos, static_cast<TDel>(varSet[j].variantSize), coverageSet[j]);
                    break;
            }
        }
    }
};

template <typename TPos, typename TSize, typename TInsBuffer>
struct VariantInfo
{
    typedef typename DeltaType::TValue TDeltaType;
    TPos refPos;
    TSize variantSize;
    TDeltaType variantType;
    TInsBuffer seqBuffer;

    VariantInfo() : refPos(0), variantSize(0), variantType(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _resolveConflicts()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline void _resolveConflicts(DeltaMap<TValue, TAlphabet, TSpec> & varStore)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TStoreIter;
    typedef typename DeltaCoverage<TDeltaMap>::Type TBitVector;
    typedef typename Position<TDeltaMap>::Type TPosition;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type TIndel;

    TStoreIter itBegin = begin(varStore, Standard());
    TStoreIter it = itBegin;
    TStoreIter itEnd = end(varStore, Standard());
    for (;it != itEnd; ++it)
    {
        if (deltaType(it) == DeltaType::DELTA_TYPE_DEL)
        {  // Resolve all variants that fall into a previously deleted region.
            TPosition endPoint = deltaDel(it) + *it;
            TStoreIter itLocal = it + 1;
            while (*itLocal < endPoint)
            {
                transform(deltaCoverage(itLocal), deltaCoverage(itLocal), deltaCoverage(it),
                          FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
                ++itLocal;
            }
        }
        if (deltaType(it) == DeltaType::DELTA_TYPE_INS)
        {  // Resolve all insertions that occur immediately before an replacement or deletion.
            TStoreIter itLocal = it + 1;
            while (*itLocal == *it)
            {
//                TMappedDelta deltaInfoInner = mappedDelta(varStore, itLocal - itBegin);
                if (deltaType(itLocal) == DeltaType::DELTA_TYPE_DEL || deltaType(itLocal) == DeltaType::DELTA_TYPE_SNP)
                {
                    TBitVector tmpVec;
                    transform(tmpVec, deltaCoverage(it), deltaCoverage(itLocal), FunctorBitwiseAnd());
                    if (!testAllZeros(tmpVec))  // At least one sequence contains ambiguous variant.
                    {
                        // Update the coverage of the nodes.
                        transform(deltaCoverage(it), deltaCoverage(it), tmpVec,
                                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
                        transform(deltaCoverage(itLocal), deltaCoverage(itLocal), tmpVec,
                                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

                        if (deltaType(itLocal) == DeltaType::DELTA_TYPE_DEL)
                        {
                            appendValue(varStore._deltaStore._indelData, TIndel(deltaDel(itLocal), deltaIns(it)));
                        }
                        else
                        {
                            TIns tmp = deltaIns(it);
                            append(tmp, deltaSnp(itLocal));
                            appendValue(varStore._deltaStore._indelData, TIndel(1, tmp));
                        }
                        // Insert the new coverage and variant info into the variant store.
                        TPosition currPos = it - begin(varStore, Standard());
                        insertValue(varStore._deltaCoverageStore._coverageData, currPos, tmpVec);
                        insertValue(varStore._deltaStore._varDataMap, currPos,
                                    (length(varStore._deltaStore._indelData) -1) | DeltaType::DELTA_TYPE_INDEL);

                        // Insert the ref position and synchronize the iterator.
                        insertValue(varStore._keys, currPos, *it);
                        it = begin(varStore, Standard()) + currPos + 1;
                        itEnd = end(varStore, Standard());  // Update end pointer.
                        itLocal = it +1;
                    }
                }
                ++itLocal;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function _extracVariants()
// ----------------------------------------------------------------------------

template <typename TVariantInfo, typename TVcfRecord>
inline int _extractVariants(String<TVariantInfo> & varString,
                            TVcfRecord & vcfRecord)
{

    StringSet<CharString> altVarSet;
    splitString(altVarSet, vcfRecord.alt, ',');

    resize(varString, length(altVarSet), Exact());

    for (unsigned i = 0; i < length(varString); ++i)  // Iterate over the variants.
    {
        CharString alt = altVarSet[i];

        int leftLcp = lcpLength(vcfRecord.ref, alt);
        CharString altR(suffix(alt, leftLcp));  // Only take from left lcp different value.
        CharString refR(suffix(vcfRecord.ref, leftLcp));  // Only take from right lcp different value.
        reverse(altR);
        reverse(refR);
        int rightLcp = lcpLength(altR, refR);

        // Extract the corresponding changes.
        Segment<CharString, InfixSegment> refInf = infix(vcfRecord.ref, leftLcp, length(vcfRecord.ref) - rightLcp);
        Segment<CharString, InfixSegment> altInf = infix(alt, leftLcp, length(alt) - rightLcp);

        TVariantInfo & info = varString[i];
        info.refPos = vcfRecord.beginPos + leftLcp;

        // Variant is a SNP or a single base deletion.
        if (length(refInf) == length(altInf))
        {
            if (length(alt) == 1 && alt[0] == '.')  // Single base deletion
            {
    //            std::cerr << "DEL ref: " << refInf << " alt: " << altInf << std::endl;
                info.variantType = DeltaType::DELTA_TYPE_DEL;
                info.variantSize = length(refInf);
            }
            else
            {
    //            std::cerr << "SNP ref: " << refInf << " alt: " << altInf << std::endl;
                info.variantType = DeltaType::DELTA_TYPE_SNP;
                info.variantSize = length(refInf);
                info.seqBuffer = altInf;
            }
        }
        // Variant is a deletion.
        else if (length(altInf) == 0)
        {
    //        std::cerr << "DEL ref: " << refInf << " alt: " << altInf << std::endl;
            info.variantType = DeltaType::DELTA_TYPE_DEL;
            info.variantSize = length(refInf);
        }
        // Variant is an Insertion
        else if (length(refInf) == 0)
        {
    //        std::cerr << "INS ref: " << refInf << " alt: " << altInf << std::endl;
            info.variantType = DeltaType::DELTA_TYPE_INS;
            info.variantSize = length(altInf);
            info.seqBuffer = altInf;
        }
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function _extractHaplotype()
// ----------------------------------------------------------------------------

template <typename TSource>
inline unsigned
_extractHaplotype(TSource const & source)
{
    // Fast switch for diploides.
    switch(source)
    {
        case '0': return 0;
        case '1': return 1;
        case '2': return 2;
        default:
            unsigned res;
            std::stringstream buffer;
            buffer << source;
            buffer >> res;
            return res;
    }
}

// ----------------------------------------------------------------------------
// Function _readVcfRecords
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet1, typename TStream>
inline int
_readVcfRecords(VcfRecordTranslator<DeltaMap<TValue, TAlphabet1> > & delegate,
                VcfIOContext & vcfContext,
                RecordReader<TStream, SinglePass<> > & vcfReader,
                JSeqHeader const & jseqHeader,
                ConverterOptions const & options)
{
    typedef DeltaMap<TValue, TAlphabet1> TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type TInsBuffer;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;
    typedef VariantInfo<TDel, TDel, TInsBuffer> TVariantInfo;
    typedef typename Value<TInsBuffer>::Type TAlphabet;


//#ifdef PRINT_DOTS
//    if (options.verbosity == 3)
    std::cerr << "\n# variants processed: ";
    static const unsigned DOT_SIZE = 10000;
    static const unsigned CHUNK_SIZE = 100000;
    unsigned pos = 0;
//#endif //PRINT_DOTS

    unsigned numSeq = length(jseqHeader._nameStore);  // Total number of sequences.
    if (options.includeReference)
        --numSeq;

    String<bool, Packed<> > altCoverage;
    resize(altCoverage, length(jseqHeader._nameStore), false, Exact());   // Coverage of all sequences.

    String<TVariantInfo> tmpInfo;
    String<String<bool, Packed<Alloc<> > > > tmpCov;
    TVariantInfo _info;
    _info.refPos = 0;
    appendValue(_info.seqBuffer, TAlphabet());
    _info.variantSize = 1;
    _info.variantType = DeltaType::DELTA_TYPE_SNP;
    appendValue(tmpInfo, _info);
    appendValue(tmpCov, altCoverage);
    // Add dummy begin node.
    delegate(tmpInfo, tmpCov, true);

    if (numSeq == 0)
        return 1;

    unsigned counter = 0;
    String<double> timeTable;
    resize(timeTable, 5, 0.0, Exact());

    CharString svFlag = "SVTYPE";

    VcfRecord vcfRecord;
    while(!atEnd(vcfReader))
    {
        ++counter;
        double timeRead = sysTime();
        if (readRecord(vcfRecord, vcfReader, vcfContext, Vcf()) != 0)
        {
            std::cerr << "Error when reading vcf file!" << std::endl;
            return 1;
        }
        timeTable[3] += sysTime() - timeRead;

//#ifdef PRINT_DOTS
        if (pos % CHUNK_SIZE == 0)
            std::cerr << " " << pos << " ";
//#endif //PRINT_DOTS

        // Only load records which are passed and have genotype information.
        toUpper(vcfRecord.filter);
        toUpper(vcfRecord.format);
        toUpper(vcfRecord.info);

        // Suppress SVs if enabled.
        if (options.suppressSVs)
        {

            if (std::search(begin(vcfRecord.format), end(vcfRecord.format), begin(svFlag), end(svFlag)) ==
                end(vcfRecord.format, Standard()))
                continue;
        }

        // Check general filter options.
        if (vcfRecord.filter == "PASS" && startsWith(vcfRecord.format, "GT") &&
            !startsWith(vcfRecord.info, "IMPRECISE") && !(vcfRecord.qual != vcfRecord.qual))
        {

            // Extracts all alternative allels for the current loci.
            String<TVariantInfo> variantInfoString;

            double timeExtract = sysTime();
            _extractVariants(variantInfoString, vcfRecord);
            timeTable[0] += sysTime() - timeExtract;

            // Coverages per alternative allel for the current loci.
            String<String<bool, Packed<> > > altCoverageSet;
            resize(altCoverageSet, length(variantInfoString), altCoverage, Exact());

            // Read the combined genotype.
            if (options.readGenotype)
            {
                for (unsigned seqId = 0; seqId < numSeq; ++seqId)
                {
                    // Process all haplotypes to generate a single genotype representing one value or not.
                    int altPos = 0;
                    for (unsigned htId = 0; htId < 2u; ++htId)  // NOTE(rmaerker): We only assume diploid.
                    {
                        double timeHaplotype = sysTime();
                        int tmpAlt = _extractHaplotype(vcfRecord.genotypeInfos[seqId][htId * 2]);  // Extract the alt position for the current sequence.
                        timeTable[1] += sysTime() - timeHaplotype;
                        if(altPos == 0 && tmpAlt != 0)
                        {
                            altPos = tmpAlt;
                            double timeAssignVec = sysTime();
                            assignValue(altCoverageSet[tmpAlt - 1], seqId, true);
                            timeTable[4] += sysTime() - timeAssignVec;
                        }
                        if (tmpAlt != 0 && tmpAlt != altPos)
                            std::cerr << "WARNING: Different allels at the same loci in: " << vcfRecord.rID  << " at "  << vcfRecord.beginPos << " for " << seqId << "!" << std::endl;
                    }
                }
            }
            else
            {
                unsigned seqId = 0;
                unsigned individualId = 0;
                while (seqId < numSeq)
                {
                    // Iterate each haplotype.
                    for (unsigned htId = 0; htId < length(options.haplotypes); ++htId, ++seqId)
                    {
                        // Extract the variant per individual per haplotype.
                        double timeHaplotype = sysTime();
                        int altPos = _extractHaplotype(vcfRecord.genotypeInfos[individualId][options.haplotypes[htId] * 2]);  // Extract the alt position for the current sequence.
                        timeTable[1] += sysTime() - timeHaplotype;

                        if (altPos != 0)
                        {
                            double timeAssignVec = sysTime();
                            assignValue(altCoverageSet[--altPos], seqId, true);
                            timeTable[4] += sysTime() - timeAssignVec;
                        }
                    }
                    ++individualId;
                }
            }

            double timeDelegate = sysTime();
            delegate(variantInfoString, altCoverageSet);
            timeTable[2] += sysTime() - timeDelegate;
            // Iterates over the variants at this position.
        }
//#ifdef PRINT_DOTS
        if (++pos % DOT_SIZE == 0)
            std::cerr << ".";
//#endif //PRINT_DOTS
    }

    if (options.verbosity >= 2)
    {
        std::cerr << "\nNumber of records: " << counter << std::endl;
        std::cerr << "Time extract variant: " << timeTable[0] << " s." << std::endl;
        std::cerr << "Time extract haplotype: " << timeTable[1] << " s." << std::endl;
        std::cerr << "Time assign Vec: " << timeTable[4] << " s." << std::endl;
        std::cerr << "Time store data: " << timeTable[2] << " s." << std::endl;
        std::cerr << "Time read data: " << timeTable[3] << " s." << std::endl;
    }
    return 0;
}

template <typename TStream, typename TRefAlphabet, typename TVarAlphabet>
int _writeToBinaryJSeqFormat(JSeqHeader & jseqHeader,
                             VcfIOContext & vcfContext,
                             RecordReader<TStream, SinglePass<> > & vcfReader,
                             ConverterOptions const & converterOptions,
                             TRefAlphabet const & /*refAlphabet*/,
                             TVarAlphabet const & /*snpAlphabet*/)
{
    typedef DeltaMap<unsigned, TVarAlphabet> TDeltaMap;

    TDeltaMap deltaMap;

    VcfRecordTranslator<TDeltaMap> delegate(deltaMap);

    double start = sysTime();
    _readVcfRecords(delegate, vcfContext, vcfReader, jseqHeader, converterOptions);
    if (converterOptions.verbosity > 1)
        std::cout << "Time for reading vcf: " << sysTime() - start << std::endl;

    _resolveConflicts(deltaMap);

    std::ofstream outputStream;
    outputStream.open(toCString(converterOptions.outputFile), std::ios_base::out);
    if (!outputStream.good())
    {
        std::cerr << "Cannot open file <"<< converterOptions.outputFile << "to write!";
        return JSeqTools::FILE_READ_ERROR;
    }

    start = sysTime();

    write(outputStream, deltaMap, jseqHeader, JSeq());
    if (converterOptions.verbosity > 1)
    {
        std::cout << "Time for writing: " << sysTime() - start << std::endl;
        std::cout << "Number of converted nodes: " << length(deltaMap);
    }

    outputStream.close();
    return 0;
}

// ----------------------------------------------------------------------------
// Function adaptVcfToJournal()
// ----------------------------------------------------------------------------

template <typename TRefAlphabet, typename TVarAlphabet>
int adaptVcfToJournalData(ConverterOptions const & converterOptions,
                          TRefAlphabet const & /*refAlphabet*/,
                          TVarAlphabet const & /*varAlphabet*/)
{
    typedef String<TRefAlphabet, Alloc<> > TReference;

    if (converterOptions.verbosity > 0)
            std::cout << "Start extracting variants.";

    // Open and read the sequence database.
    std::ifstream vcfStream;
    vcfStream.open(toCString(converterOptions.inputFile), std::ios_base::in);
    if (!vcfStream.good())
    {
        std::cerr << "ERROR: Could not open "<<  converterOptions.inputFile << "!" << std::endl;
        return JSeqTools::FILE_READ_ERROR;
    }
    VcfHeader vcfHeader;
    VcfIOContext vcfContext(vcfHeader.sequenceNames, vcfHeader.sampleNames);
    RecordReader<std::ifstream, SinglePass<> > vcfReader(vcfStream);
    read(vcfHeader, vcfReader, vcfContext, Vcf());

    unsigned numSamples = _min(length(vcfHeader.sampleNames), converterOptions.numIndividuals);
    // Resize to the maximum number of individuals that can be parsed.
    unsigned totalSequences = numSamples * length(converterOptions.haplotypes);

    JSeqHeader jseqHeader;

    if (IsSameType<TVarAlphabet, Dna>::VALUE || IsSameType<TVarAlphabet, Rna>::VALUE)
        jseqHeader._fileInfos._snpCompression = true;
    else
        jseqHeader._fileInfos._snpCompression = false;

    // Prepare the sequence names for the output of the delta file.
    if (converterOptions.readGenotype)
    {
        resize(jseqHeader._nameStore, numSamples, Exact());
        for (unsigned i = 0; i < numSamples; ++i)
            assignValue(jseqHeader._nameStore, i, getValue(vcfHeader.sampleNames, i));
    }
    else
    {
        resize(jseqHeader._nameStore, totalSequences, Exact());
        unsigned seqId = 0;
        for (unsigned individualId = 0; individualId < numSamples; ++individualId)
            for (unsigned j = 0; j < length(converterOptions.haplotypes); ++j, ++seqId)
            {
                CharString fileName = getValue(vcfHeader.sampleNames, individualId);
                std::stringstream nameStream;
                nameStream << "_ht_" << converterOptions.haplotypes[j];
                append(fileName, nameStream.str());
                assignValue(jseqHeader._nameStore, seqId, fileName);
            }
    }

    if (converterOptions.includeReference)
        appendValue(jseqHeader._nameStore, "reference");

    // Load the reference file.
    TReference reference;
    if (_loadSequenceFasta(jseqHeader._refInfos._refId, reference, converterOptions.vcfReferenceFile) != JSeqTools::FILE_READ_OK)
        return JSeqTools::FILE_READ_ERROR;

    jseqHeader._refInfos._refFile = converterOptions.vcfReferenceFile;

    _writeToBinaryJSeqFormat(jseqHeader, vcfContext, vcfReader, converterOptions, TRefAlphabet(), TVarAlphabet());


    vcfStream.close();
    if (converterOptions.verbosity > 0)
        std::cout << "\tDone!" << std::endl;

    return 0;
}

}

#endif // EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_ADAPT_VCF_H_
