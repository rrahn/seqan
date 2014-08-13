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
// Facade header for the jseq tools.
// ==========================================================================

#ifndef EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_FACADE_HEADER_H_
#define EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_FACADE_HEADER_H_

#include <iostream>
#include <algorithm>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/journaled_set.h>

#include <seqan/basic/basic_view.h>  // TODO(rmaerker): Probably no global header.
#include <seqan/find.h>  // TODO(rmaerker): Probably no global header.
#include <seqan/parallel.h>

#include <seqan/seq_io.h>
#include <seqan/journaled_string_tree.h>
#include <seqan/find_journaled_string_tree.h>

//#include "../gdf_tools/gdf_io_base.h"
//#include "../gdf_tools/gdf_io_header.h"
//#include "../gdf_tools/gdf_io_read.h"
//#include "../gdf_tools/gdf_io_write.h"

#include "jst_bench_options.h"
#include "jst_bench_base.h"
#include "jst_bench_io.h"

#endif // EXTRAS_APPS_JSEQ_TOOLS_JSEQ_TOOLS_FACADE_HEADER_H_
