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
// Tool to convert alignment formats into gdf format.
// ==========================================================================

#define PRINT_DOTS

#include "gdf_tools_facade_header.h"
#include "gdf_converter_vcf.h"  //  Converting vcf files into gdf format.

using namespace seqan;

// ----------------------------------------------------------------------------
// Function readFinderOptions()
// ----------------------------------------------------------------------------

inline ArgumentParser::ParseResult
readAdaptOptions(ConverterOptions & options, ArgumentParser const & parser)
{
    // Read the main options of the program.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    // Read gdf file from argument list.
    getArgumentValue(options.outputFile, parser, 0);
//    getOptionValue(options.outputFile, parser, "gdf-file");

    // Read VCF Options.
    getOptionValue(options.inputFile, parser, "vcf-file");
    getOptionValue(options.vcfReferenceFile, parser, "vcf-reference");

    CharString tmpAlphabetRef;
    getOptionValue(tmpAlphabetRef, parser, "ra");
    if (tmpAlphabetRef == "dna")
        options.refAlphabet = 0;
    if (tmpAlphabetRef == "dna5")
        options.refAlphabet = 1;

    CharString tmpAlphabetVar;
    getOptionValue(tmpAlphabetVar, parser, "va");
    if (tmpAlphabetVar == "dna")
        options.varAlphabet = 0;
    if (tmpAlphabetVar == "dna5")
        options.varAlphabet = 1;

    if (isSet(parser, "n"))
        getOptionValue(options.numIndividuals, parser, "n");
    else
        options.numIndividuals = - 1;

    if (isSet(parser, "genotype"))
        options.readGenotype = true;
    else  // Read haplotype information.
    {
        std::stringstream tmpStr;
        String<std::string> tmpOptions = getOptionValues(parser, "haplotype");
        resize(options.haplotypes, length(tmpOptions), Exact());

        for (unsigned i = 0; i < length(tmpOptions); ++i)
        {
            tmpStr.str("");
            tmpStr.clear();
            tmpStr << tmpOptions[i];
            tmpStr >> options.haplotypes[i];
        }
        std::sort(begin(options.haplotypes, Standard()), end(options.haplotypes, Standard()));
    }

//    if (isSet(parser, "include-reference"))
//        options.includeReference = true;

    if(isSet(parser, "suppress-sv"))
        options.suppressSVs = true;

//    for (unsigned j = 0; j < length(options.haplotypes); ++j)
//        std::cerr << options.haplotypes[j] << "\t";
//    std::cerr << std::endl;

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function _setupAndParseArguments()
// ----------------------------------------------------------------------------

template <typename TArgCount, typename TArgValues>
ArgumentParser::ParseResult
_setupAndParseArguments(ConverterOptions & options,
                        TArgCount argc,
                        TArgValues const & argv)
{
    ArgumentParser parser("gdf_converter");
    // Set short description, version, and date.
    setShortDescription(parser, "GDF-Converter for sequence alignments.");
    setVersion(parser, "0.1");
    setDate(parser, "April 2014");

    // Define usage line and long description.
    addUsageLine(parser," \"[\\fIOPTIONS\\fP] \\fIGDF-FILE\\fP\" ");
    addDescription(parser, "The GDF-Converter converts alignment formats into the compressed genome delta format (GDF)."
                           "The GDF format solely stores the differences of all sequences to a common "
                           "reference sequence and a bit vector for each variant to be mask the sequences "
                           "that share this variant.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE, "GDF-FILE"));

    addSection(parser, "General Options");
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

//    addOption(parser, seqan::ArgParseOption("g", "gdf-file", "Name of the GDF file.", ArgParseArgument::OUTPUTFILE));
//    setValidValues(parser, "gdf-file", "gdf");

    addSection(parser, "VCF Conversion");
    addOption(parser, ArgParseOption("vcf", "vcf-file", "Vcf file to be converted.", ArgParseArgument::INPUTFILE));
    setValidValues(parser, "vcf-file", "vcf");

    addOption(parser, seqan::ArgParseOption("r", "vcf-reference", "Fasta file of the reference sequences.", ArgParseArgument::INPUTFILE));
    setValidValues(parser, "vcf-reference", "fa fasta");

//    addOption(parser, ArgParseOption("", "include-reference", "Includes the reference within the gdf file."));

    addOption(parser, ArgParseOption("ra", "reference-alphabet", "Alphabet of the reference sequence.", ArgParseArgument::STRING));
    setValidValues(parser, "ra", "dna dna5");
    setDefaultValue(parser, "ra", "dna5");

    addOption(parser, ArgParseOption("va", "vcf-alphabet", "Alphabet of the variants.", ArgParseArgument::STRING));
    setValidValues(parser, "va", "dna dna5");
    setDefaultValue(parser, "va", "dna");

    addOption(parser, ArgParseOption("ht", "haplotype", "Treats each haplotype as a single sequence. Repeat this option to specify multiple haplotypes.", ArgParseArgument::INTEGER, "NUM", true));
    setDefaultValue(parser, "ht", "0");

    addOption(parser, ArgParseOption("n", "firstn", "If specified only the first \"n\" sequences are converted. Takes all individuals per default", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "n", "1");

    addOption(parser, ArgParseOption("", "genotype", "Generates genotype by merging the haplotypes into a single sequence. Note, this option will overwrite the option \"haplotype\" if set."));

    addOption(parser, ArgParseOption("", "suppress-sv", "Suppresses conversion of structural variants."));

    // Add Examples Section.
//    addTextSection(parser, "Examples");
//    addListItem(parser, "\\fBcolossus\\fP \\fB-v\\fP \\fItext\\fP",
//                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        return res;
    }


    return readAdaptOptions(options, parser);
}

template <typename TRefAlphabet>
inline int adaptVcfToJournalData(ConverterOptions const & options, TRefAlphabet /*alph*/)
{
    switch (options.varAlphabet)
    {
        case 0: return adaptVcfToJournalData(options, TRefAlphabet(), Dna());
        case 1: return adaptVcfToJournalData(options, TRefAlphabet(), Dna5());
        default: return -1;
    }
}

inline int adaptVcfToJournalData(ConverterOptions const & options)
{
    switch (options.refAlphabet)
    {
        case 0: return adaptVcfToJournalData(options, Dna());
        case 1: return adaptVcfToJournalData(options, Dna5());
        default: return -1;
    }
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const** argv)
{

    ConverterOptions options;
    ArgumentParser::ParseResult res = _setupAndParseArguments(options, argc, argv);

    if (res == ArgumentParser::PARSE_HELP || res == ArgumentParser::PARSE_VERSION)
        return res;

    if (res != ArgumentParser::PARSE_OK)
    {
        std::cerr << "Errors while parsing the input occurred!" << std::endl;
        return res;
    }

    // We need to check the file type.
    if (endsWith(options.inputFile, "vcf") || endsWith(options.inputFile, "VCF"))
    {
        return adaptVcfToJournalData(options);
    }
    std::cerr << "Unsupported alignment file!" << std::endl;
    return 0;
}
