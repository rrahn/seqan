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
// Implements the pattern state for the horspool algorithm.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_IMPL_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_IMPL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TContainer, typename TPattern, typename TSpec>
struct Finder2<TContainer, TPattern, DataParallel<TSpec> >
{
    typedef typename FinderFunctor<Finder2>::Type TFinderFunctor;

    TContainer*     _containerPtr;
    TFinderFunctor  _finderFunctor;

    bool _needReinit;

    Finder2() : _containerPtr(NULL), _needReinit(true)
    {}


    Finder2(TContainer & hystk) : _containerPtr(&hystk), _finderFunctor(), _needReinit(true)
    {}

    // Copy constructor.
    Finder2(Finder2 const & other) : _containerPtr(other._containerPtr),
                                     _finderFunctor(other._finderFunctor),
                                     _needReinit(other._needReinit)
    {}

    Finder2 & operator=(Finder2 const & other)
    {
        if (this != &other)
        {
            _containerPtr = other._containerPtr;
            _finderFunctor = other._finderFunctor;
            _needReinit = other._needReinit;
        }
        return *this;
    }

//    Finder2(TContainer & hystk, TPattern & pattern) :
//            _containerPtr(&hystk),
//            _finderFunctor(pattern),
//            _needReinit(false)
//    {}
//
//    template <typename TMinScore>
//    Finder2(TContainer & hystk, TPattern & pattern, TMinScore const & minScore) :
//            _containerPtr(&hystk),
//            _finderFunctor(pattern, minScore),
//            _needReinit(false)
//    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TContainer, typename TPattern, typename TSpec>
struct GetTraverserForFinder_<Finder2<TContainer, TPattern, TSpec> const>
{
    typedef typename PatternSpecificTraversalSpec<TPattern>::Type TTraversalSpec_;
    typedef Traverser<TContainer const, TTraversalSpec_> Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct GetTraverserForFinder_<Finder2<TContainer, TPattern, TSpec> >
{
    typedef typename PatternSpecificTraversalSpec<TPattern>::Type TTraversalSpec_;
    typedef Traverser<TContainer, TTraversalSpec_> Type;
};


// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct Position<Finder2<TContainer, TPattern, DataParallel<TSpec> > >
{
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder_;
    typedef typename GetTraverserForFinder_<TFinder_>::Type TTraverser_;

    typedef typename Position<TTraverser_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct Position<Finder2<TContainer, TPattern, DataParallel<TSpec> > const>
{
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder_;
    typedef typename GetTraverserForFinder_<TFinder_>::Type TTraverser_;

    typedef typename Position<TTraverser_ const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct Container<Finder2<TContainer, TPattern, DataParallel<TSpec> > >
{
    typedef TContainer Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct Container<Finder2<TContainer, TPattern, DataParallel<TSpec> > const>
{
    typedef TContainer const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function compute()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TIterator>
inline typename ComputeState<TContainer>::Type
compute(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder, TIterator const & iter)
{
    typedef typename PatternSpecificTraversalSpec<TPattern>::Type TTraversalSpec;
    typedef Traverser<TContainer, TTraversalSpec> TTraverser;

    typedef typename ComputeState<TTraverser>::Type TComputeState;

    TComputeState state(false, 1);
    finder._finderFunctor(state, iter);
    return state;
}

// ----------------------------------------------------------------------------
// Function getPattern()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TIterator>
inline TPattern &
getPattern(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder)
{
    return finder._finderFunctor._pattern;
}

template <typename TContainer, typename TPattern, typename TSpec, typename TIterator>
inline TPattern const &
getPattern(Finder2<TContainer, TPattern, DataParallel<TSpec> > const & finder)
{
    return finder._finderFunctor._pattern;
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
inline typename Container<Finder2<TContainer, TPattern, DataParallel<TSpec> > >::Type &
container(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder)
{
    return *finder._containerPtr;
}

template <typename TContainer, typename TPattern, typename TSpec>
inline typename Container<Finder2<TContainer, TPattern, DataParallel<TSpec> > const>::Type &
container(Finder2<TContainer, TPattern, DataParallel<TSpec> > const & finder)
{
    return *finder._containerPtr;
}


template <typename TContainer, typename TPattern, typename TSpec, typename TScoreLimit>
inline void
_init(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
      TPattern & pattern,
      TScoreLimit const & scoreLimit,
      True const & /*errorsSupported*/)
{
    _init(finder._finderFunctor, pattern, scoreLimit);
}

template <typename TContainer, typename TPattern, typename TSpec, typename TScoreLimit>
inline void
_init(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
      TPattern & pattern,
      TScoreLimit const & /*scoreLimit*/,
      False const & /*errorsSupported*/)
{
    _init(finder._finderFunctor, pattern);
}

// ----------------------------------------------------------------------------
// Function find()                                              [Parallel case]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate, typename TParallelTag>
inline void
find(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
     TPattern & pattern,
     TDelegate & delegate,
     int scoreLimit = 0,
     Tag<TParallelTag> tag = Serial())    // A tag must be specified.
{
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder;
    typedef typename FinderFunctor<TFinder>::Type TFinderFunctor;
    typedef typename GetTraverserForFinder_<TFinder>::Type TTraverser;
    typedef typename VariantData<TContainer>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap, Rooted>::Type TDeltaIter;
    typedef typename Position<TContainer>::Type TPosition;


    requireJournal(container(finder), tag);  // Build the journal set on demand - works in parallel too.

    if (finder._needReinit)
    {
        // finder._finderFunctor(pattern, initState);
        _init(finder, pattern, scoreLimit, typename ErrorsSupported<TFinderFunctor>::Type());
        finder._needReinit = false;
    }

    // Initialize the traverser.
    TTraverser traverser(container(finder), length(needle(pattern)) - scoreLimit);
    // Parallelize over set of branch nodes. -> Maybe not most efficient.
    Splitter<TDeltaIter> nodeSplitter(begin(variantData(container(finder)), Rooted()), end(variantData(container(finder)), Rooted()), Parallel());

    StringSet<String<TPosition> > mergePointOverlaps;
    resize(mergePointOverlaps, length(nodeSplitter));

    SEQAN_OMP_PRAGMA(parallel for firstprivate(traverser))
    for (unsigned jobId = 0; jobId < length(nodeSplitter); ++jobId)
    {
        TFinder threadFinder = finder;  // Copy basic initialized finder.  // NOTE(rmaerker): Copying the states is probably much faster than do a initialization for each finder per thread.

//        printf("In thread: %i of %i\n", omp_get_thread_num(), omp_get_num_threads());
        TPosition hostBeginPos = 0;
        if (jobId != 0u)
            hostBeginPos = *nodeSplitter[jobId] - (windowSize(traverser) - 1);
        TPosition hostEndPos = *nodeSplitter[jobId + 1];
        if (jobId == static_cast<unsigned>(omp_get_num_threads() - 1))
            hostEndPos = length(host(container(finder)));

        _initSegment(traverser, nodeSplitter[jobId], nodeSplitter[jobId + 1], hostBeginPos, hostEndPos);
        traverse(traverser, threadFinder, delegate);

        // TODO(rmaerker): Smooth the results in case we found a hit that was deleted through a deletion of a previous delta.
        if (length(traverser._mergePointStack._mergePoints) > 1u)
            mergePointOverlaps[jobId] = traverser._mergePointStack._mergePoints;
    }

    // We might need to update some hits here.
    // How can we return this?
}

// ----------------------------------------------------------------------------
// Function find()                                                [Serial case]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
     TPattern & pattern,
     TDelegate & delegate,
     int scoreLimit = 0)                // Use default all threads available.
{
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder;
    typedef typename FinderFunctor<TFinder>::Type TFinderFunctor;
    typedef typename GetTraverserForFinder_<TFinder>::Type TTraverser;

    if (finder._needReinit)
    {
        // finder._finderFunctor(pattern, initState);
        _init(finder, pattern, scoreLimit, typename ErrorsSupported<TFinderFunctor>::Type());
        finder._needReinit = false;
    }

    requireJournal(container(finder), Serial());

    // Set up the string tree traversal.
    TTraverser traverser(container(finder), length(needle(pattern)) - scoreLimit);
    traverse(traverser, finder, delegate);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_DATA_PARALLEL_FIND_DATA_PARALLEL_IMPL_H_
