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
//  Implements the new interface for calling alingment algorithms.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_

namespace seqan {
namespace impl {
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct SimdWorker
{
    template <typename TQueueContext>
    inline void
    operator()(TQueueContext & queueContext)
    {
        using TTask = decltype(popFront(queueContext.mQueue));

        lockWriting(queueContext.mQueue);
        std::vector<TTask> tasks;
        while (true)
        {
            TTask task = nullptr;
            tasks.clear();
            {
                std::lock_guard<decltype(queueContext.mLock)> scopedLock(queueContext.mLock);
                task = popFront(queueContext.mQueue);
                if (task == nullptr)
                    return;

                if (length(queueContext.mQueue) >= TQueueContext::VECTOR_SIZE - 1)
                {
                    for (unsigned i = 0; i < TQueueContext::VECTOR_SIZE - 1; ++i)
                        tasks.push_back(popFront(*workQueuePtr));
                }
            }

            SEQAN_ASSERT(task != nullptr);
            task->template execute(*workQueuePtr, tasks, mThreadId);
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Now we can implement different strategies to compute the alignment.
template <typename TScore, typename TDPTraits, typename TExecutionTraits,
          typename ...Ts>
void align_batch(DPConfig<TScore, TDPTraits, TExecutionTraits> const & config,
                 Ts ...&& args)
{
    BatchAlignmentExecutor<typename TExecutionTraits::TParallelPolixy, typename TExecutionTraits::TSchedulingPolicy>::run(config, std::forward<Ts>(args)...);
}

}  // namespace impl

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_
