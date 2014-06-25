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
// Implements hamming distance based comparison of two strings.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSize>
struct StateHammingSimple_
{
    TSize errorCount;
    TSize maxErrors;
};

template <typename TPattern>
class FinderExtensionPoint<TPattern,  HammingSimple>
{
public:
    typedef typename GetState<FinderExtensionPoint>::Type TState;

    TState state;

    FinderExtensionPoint()
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition                               [Hamming]
// ----------------------------------------------------------------------------

template <typename T>
struct ContextIteratorPosition<FinderExtensionPoint<T, HammingSimple> >
{
    typedef ContextPositionRight Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext                                    [Hamming]
// ----------------------------------------------------------------------------

template <typename T>
struct RequireFullContext<FinderExtensionPoint<T, HammingSimple> > : False{};

// ----------------------------------------------------------------------------
// Metafunction RegisteredExtensionPoint                              [Hamming]
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct RegisteredExtensionPoint<Pattern<TNeedle, HammingSimple> >
{
    typedef Pattern<TNeedle, HammingSimple> TPattern_;
    typedef FinderExtensionPoint<TPattern_, HammingSimple> Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState                                              [Hamming]
// ----------------------------------------------------------------------------

template <typename T>
struct GetState<FinderExtensionPoint<T, HammingSimple> >
{
    typedef StateHammingSimple_ Type;
};

template <typename T>
struct GetState<FinderExtensionPoint<T, HammingSimple> const>
{
    typedef StateHammingSimple_ const Type;
};

// ============================================================================
// Functions
// ============================================================================

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_EXTENSION_HAMMING_H_
