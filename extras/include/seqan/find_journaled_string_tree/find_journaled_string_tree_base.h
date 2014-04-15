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
// Tags and structures used globally for journaled string tree finder.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_BASE_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag BitAlgorithmLongNeedle
// ----------------------------------------------------------------------------

struct BitAlgorithmLongNeedle_;
typedef Tag<BitAlgorithmLongNeedle_> BitAlgorithmLongNeedle;

// ----------------------------------------------------------------------------
// Tag BitAlgorithmSmallNeedle
// ----------------------------------------------------------------------------

struct BitAlgorithmSmallNeedle_;
typedef Tag<BitAlgorithmSmallNeedle_> BitAlgorithmSmallNeedle;

// ----------------------------------------------------------------------------
// Tag DataParallel
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct DataParallel;

// ----------------------------------------------------------------------------
// Struct FinderFunctor
// ----------------------------------------------------------------------------

template <typename TFinder>
struct FinderFunctor{};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition
// ----------------------------------------------------------------------------

template <typename T>
struct ContextIteratorPosition
{
    typedef ContextPositionLeft Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext
// ----------------------------------------------------------------------------

template <typename T>
struct RequireFullContext : True{};

// ----------------------------------------------------------------------------
// Metafunction TraversalSpec
// ----------------------------------------------------------------------------

template <typename T>
struct GetTraverserForFinder
{
    typedef TraverserSpec<> Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetTraverserForFinder_
// ----------------------------------------------------------------------------

template <typename TFinder>
struct GetTraverserForFinder_;

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_BASE_H_
