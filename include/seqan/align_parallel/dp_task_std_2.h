// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
//  DP Task implementation.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_STD_2_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_STD_2_H_

namespace seqan
{
namespace impl
{
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename ...TArgs>
struct DPTaskStdContext
{
    TTaskQueue &        mQueue;
    TThreadLocal &      mThreadLocal;
    TTileBuffer &       mBuffer;
    TSeqBlocksH const & mSeqBlocksH;
    TSeqBlocksV const & mSeqBlocksV;
    TDPConfig const &   mConfig;
    uint16_t            mAlignInstanceId;

    // Notification.
    std::mutex          lock;
    std::conditional    readyEvent;
    bool                ready;
};

template <typename TTaskContext>
class DPTaskStd : public DPTaskBase<DPTaskStd<TTaskContext>>
{
    TTaskContext& mContext;      // A set of shared data.
    bool          mLastTile;
    // Context can keep also the simd score.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TTaskContext>
inline void
execute(DPTaskStd & task)
{
    inline void
    executeScalar(task);
    if (task.mLastTile)
    {
        {
            std::lock_guard<std::mutex> lck(task.mContext.lock);
            task.mContext.ready = true;
        }
        task.mContext.readyEvent.notify_one();
    }
}  // namespace impl
}  // namespace seqan
#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_STD_2_H_
