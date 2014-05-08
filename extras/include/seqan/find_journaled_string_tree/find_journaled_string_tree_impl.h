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
// Implements the finder interface for the journaled string tree.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_IMPL_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_IMPL_H_

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
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder;
    typedef typename FinderExtension<TFinder>::Type TExtensionFunctor;

    TContainer*         _containerPtr;
    TExtensionFunctor   _extensionFunctor;
    bool                _needReinit;

    Finder2() : _containerPtr(NULL), _needReinit(true)
    {}


    Finder2(TContainer & hystk) : _containerPtr(&hystk), _extensionFunctor(), _needReinit(true)
    {}

    // Copy constructor.
    Finder2(Finder2 const & other) : _containerPtr(other._containerPtr),
                                     _extensionFunctor(other._extensionFunctor),
                                     _needReinit(other._needReinit)
    {}

    // Assignment Operator
    Finder2 & operator=(Finder2 const & other)
    {
        if (this != &other)
        {
            _containerPtr = other._containerPtr;
            _extensionFunctor = other._extensionFunctor;
            _needReinit = other._needReinit;
        }
        return *this;
    }
};

template <typename TContainer, typename TPattern, typename TSpec>
SEQAN_CONCEPT_IMPL(Finder2<TContainer, TPattern, DataParallel<TSpec> >,       (JstTraversalConcept));

template <typename TContainer, typename TPattern, typename TSpec>
SEQAN_CONCEPT_IMPL(Finder2<TContainer, TPattern, DataParallel<TSpec> > const, (JstTraversalConcept));

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct ContextIteratorPosition<Finder2<TContainer, TPattern, TSpec> >
{
    typedef Finder2<TContainer, TPattern, TSpec> TFinder_;
    typedef typename FinderExtension<TFinder_>::Type TFinderFunctor_;
    typedef typename ContextIteratorPosition<TFinderFunctor_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct ContextIteratorPosition<Finder2<TContainer, TPattern, TSpec> const > :
    ContextIteratorPosition<Finder2<TContainer, TPattern, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct RequireFullContext<Finder2<TContainer, TPattern, TSpec> >
{
    typedef Finder2<TContainer, TPattern, TSpec> TFinder_;
    typedef typename FinderExtension<TFinder_>::Type TFinderFunctor_;
    typedef typename RequireFullContext<TFinderFunctor_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct RequireFullContext<Finder2<TContainer, TPattern, TSpec> const > :
    RequireFullContext<Finder2<TContainer, TPattern, TSpec> >{};

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

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct GetState<Finder2<TContainer, TPattern, DataParallel<TSpec> > >
{
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder_;
    typedef typename FinderExtension<TFinder_>::Type TFinderExtension_;
    typedef typename ExtensionState<TFinderExtension_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct GetState<Finder2<TContainer, TPattern, DataParallel<TSpec> > const> :
    GetState<Finder2<TContainer, TPattern, DataParallel<TSpec> > >{};

// ----------------------------------------------------------------------------
// Metafunction GetTraverserForFinder_
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct GetJstTraverser<Finder2<TContainer, TPattern, TSpec> >
{
    typedef Finder2<TContainer, TPattern, TSpec> TFinder_;
    typedef typename ContextIteratorPosition<TFinder_>::Type TContextPosition_;
    typedef typename RequireFullContext<TFinder_>::Type TRequireContext_;
    typedef typename GetState<TFinder_>::Type TState;
    typedef JstTraverser<TContainer, TState, JstTraverserSpec<TContextPosition_, TRequireContext_> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getExtensionState()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
inline typename ExtensionState<ExtensionFunctor<TFinder, TSpec> >::Type
getExtensionState(ExtensionFunctor<TFinder, TSpec> const & /*extensionFunctor*/)
{
    return typename ExtensionState<ExtensionFunctor<TFinder, TSpec> >::Type();
}

// ----------------------------------------------------------------------------
// Function setExtensionState()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec, typename TState>
inline void
setExtensionState(ExtensionFunctor<TFinder, TSpec> const & /*extensionFunctor*/,
                  TState const & /*state*/)
{
    // no-op function.
}

// ----------------------------------------------------------------------------
// Function getState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
inline typename GetState<Finder2<TContainer, TPattern, DataParallel<TSpec> > const>::Type
getState(Finder2<TContainer, TPattern, DataParallel<TSpec> > const & finder)
{
    return getExtensionState(finder._extensionFunctor);  // Delegates to extension functor.
}

// ----------------------------------------------------------------------------
// Function setState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TState>
inline void
setState(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
         TState const & state)
{
    setExtensionState(finder._extensionFunctor, state);  // Delegates to extension functor.
}

// ----------------------------------------------------------------------------
// Function execute()
// ----------------------------------------------------------------------------

template <typename TResult, typename TFinder, typename TExtensionSpec, typename TContextIter>
inline void
execute(TResult & res,
        ExtensionFunctor<TFinder, TExtensionSpec> & extensionFunctor,
        TContextIter & contextIter)
{
    extensionFunctor(res, contextIter);
}

// ----------------------------------------------------------------------------
// Function deliverContext()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate, typename TTraverser, typename TTag>
inline typename Size<TTraverser>::Type
deliverContext(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
               TDelegate & delegateFunctor,
               TTraverser & traverser,
               TTag const & /*traverserState*/)
{
    typedef typename Size<TContainer>::Type TSize;

    Pair<bool, TSize> res(false, 1);
    execute(res, finder._extensionFunctor, contextIterator(traverser, TTag()));

#ifdef DEBUG_DATA_PARALLEL
    if (res.i1)  // Interrupt: Call the DelegateFunctor.
    {
        _printContext(traverser);
        delegateFunctor(traverser);
    }
#else
    if (res.i1)  // Interrupt: Call the DelegateFunctor.
        delegateFunctor(traverser);
#endif
    // Return to the traverser and continue.
    return res.i2;
}

// ----------------------------------------------------------------------------
// Function getPattern()
// ----------------------------------------------------------------------------
// TODO(rmaerker): Where do we need this?
template <typename TContainer, typename TPattern, typename TSpec, typename TIterator>
inline TPattern &
getPattern(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder)
{
    return finder._extensionFunctor._pattern;
}

template <typename TContainer, typename TPattern, typename TSpec, typename TIterator>
inline TPattern const &
getPattern(Finder2<TContainer, TPattern, DataParallel<TSpec> > const & finder)
{
    return finder._extensionFunctor._pattern;
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

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TExtensionSpec, typename TPattern, typename TScoreLimit>
inline void
init(ExtensionFunctor<TFinder, TExtensionSpec> & extensionFunctor,
     TPattern & pattern,
     TScoreLimit const & /*scoreLimit*/)
{
    init(extensionFunctor, pattern);
}

// ----------------------------------------------------------------------------
// Function find()                                                     [Serial]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder2<TContainer, TPattern, DataParallel<TSpec> > & finder,
     TPattern & pattern,
     TDelegate & delegate,
     int scoreLimit = 0)
{
    typedef Finder2<TContainer, TPattern, DataParallel<TSpec> > TFinder;
    typedef typename GetJstTraverser<TFinder>::Type TTraverser;

    // Set up the journaled string tree traversal.
    TTraverser traverser(container(finder), length(needle(pattern)) - scoreLimit);
#ifdef PROFILE_JST_INTERN
    double timeBuild = sysTime();
#endif
    while (journalNextBlock(container(finder), contextSize(traverser)))  // Generating the journal index for the next block.
    {
        if (finder._needReinit)
        {
            init(finder._extensionFunctor, pattern, scoreLimit);
            finder._needReinit = false;
        }

        _reinitBlockEnd(traverser);
#ifdef PROFILE_JST_INTERN
        std::cerr << "Time-C: " << sysTime() - timeBuild << " s" << std::endl;
        double timeSearch = sysTime();
#endif
        traverse(finder, delegate, traverser);
#ifdef PROFILE_JST_INTERN
        std::cerr << "Time-S: " << sysTime() - timeSearch << " s" << std::endl;
        timeBuild = sysTime();
#endif
    }
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_IMPL_H_
