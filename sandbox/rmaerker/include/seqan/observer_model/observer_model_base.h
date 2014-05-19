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
// Implements a generic observable model.
// ==========================================================================

#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_OBSERVER_MODEL_OBSERVER_MODEL_BASE_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_OBSERVER_MODEL_OBSERVER_MODEL_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TObserver>
struct Observable
{
    TObserver * _observerPtr;

    Observable() : _observerPtr((TObserver *) 0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TObserver>
inline void addObserver(Observable<void> & observable, TObserver & /*obs*/)
{
    // no-op.
}

template <typename TObserver>
inline void addObserver(Observable<TObserver> & observable, TObserver & observer)
{
    observable._observerPtr = &observer;
}

// ----------------------------------------------------------------------------
// Function deleteObserver()
// ----------------------------------------------------------------------------

template <typename TObserver>
inline void
deleteObserver(Observable<void> & observable)
{
    //no-op
}

template <typename TObserver>
inline void
deleteObserver(Observable<TObserver> & observable)
{
    observable._observerPtr = (TObserver*)0;
}

// ----------------------------------------------------------------------------
// Function observer()
// ----------------------------------------------------------------------------

inline Nothing observer(Observable<void> & observable)
{
    return Nothing();
}

template <typename TObserver>
inline TObserver &
observer(Observable<TObserver> & observable)
{
    return *observable._observerPtr;
}

// ----------------------------------------------------------------------------
// Function update()
// ----------------------------------------------------------------------------

// Supports, default events for several parameter lists.
template <typename TObserver, typename TParam1>
inline void
update(Observable<TObserver> const &, TParam1 const &)
{
    // no-op.
}

template <typename TObserver, typename TParam1, typename TParam2>
inline void
update(Observable<TObserver> const &, TParam1 const &, TParam2 const &)
{
    // no-op.
}

template <typename TObserver, typename TParam1, typename TParam2, typename TParam3>
inline void
update(Observable<TObserver> const &, TParam1 const &, TParam2 const &, TParam3 const &)
{
    // no-op.
}

template <typename TObserver, typename TParam1, typename TParam2, typename TParam3, typename TParam4>
inline void
update(Observable<TObserver> const &, TParam1 const &, TParam2 const &, TParam3 const &, TParam4 const &)
{
    // no-op.
}

template <typename TObserver, typename TParam1, typename TParam2, typename TParam3, typename TParam4, typename TParam5>
inline void
update(Observable<TObserver> const &, TParam1 const &, TParam2 const &, TParam3 const &, TParam4 const &, TParam5 const &)
{
    // no-op.
}

template <typename TObserver, typename TParam1, typename TParam2, typename TParam3, typename TParam4, typename TParam5,
        typename TParam6>
inline void
update(Observable<TObserver> const &, TParam1 const &, TParam2 const &, TParam3 const &, TParam4 const &, TParam5 const &,
       TParam6 const &)
{
    // no-op.
}

#define SEQAN_NOTIFY_OBSERVER(observable, ...)  update(observer(observable), __VA_ARGS__);

}  // namespace seqan

#endif  // SANDBOX_RMAERKER_INCLUDE_SEQAN_OBSERVER_MODEL_OBSERVER_MODEL_BASE_H_
