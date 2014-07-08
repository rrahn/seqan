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
// Facade header for journaled string tree finder.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/journaled_set.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/journaled_string_tree.h>

// ===========================================================================
// Basics.
// ===========================================================================

#include <seqan/find_journaled_string_tree/find_extension_concept.h>
#include <seqan/find_journaled_string_tree/find_journaled_string_tree_base.h>
#include <seqan/find_journaled_string_tree/find_journaled_string_tree_impl.h>

// ===========================================================================
// Exact Search.
// ===========================================================================

#include <seqan/find_journaled_string_tree/find_extension_simple.h>
#include <seqan/find_journaled_string_tree/find_extension_horspool.h>
#include <seqan/find_journaled_string_tree/find_extension_shift_and.h>
#include <seqan/find_journaled_string_tree/find_extension_shift_or.h>

// ===========================================================================
// Approximate Search.
// ===========================================================================

#include <seqan/find_journaled_string_tree/find_extension_myers_ukkonen.h>


#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_
