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

// TODO(rmaerker): Finder_ class needs to be documented.
/*!
 * @class JstFinder
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @implements JstTraversalConcept
 * @brief Searches for pattern matches by traversing multiple sequences simultaneously using a
 * @link JournaledStringTree @endlink.
 *
 * The data parallel finder searches a pattern in a set of sequences simultaneously using a @link JournaledStringTree
 * @endlink as the haystack and by triggering the traversal. It implements the @link JstTraversalConcept @endlink
 * in order to evaluate the different sequence contexts explored during the traversal. The finder itself can be
 * seen as a register for extension points. By registering the @link FinderExtensionPoint @endlink using the metafunction
 * @link JstFinder#RegisteredExtensionPoint @endlink any algorithm can be plugged into the interface generically in order to
 * search multiple sequences simultaneously.
 *
 * @signature template <typename THaystack, typename TPattern>
 *            struct Finder_<THaystack, TPattern, JstFinder>;
 *
 * @tparam  THaystack   The type of the haystack to be searched for a pattern match. Of type @link JournaledStringTree @endlink.
 * @tparam  TPattern    The type of the pattern used for searching. Of type @link Pattern @endlink.
 */

/*!
 * @fn JstFinder::JstFinder
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief constructor
 *
 * @signature JstFinder();
 * @signature JstFinder(haystack);
 * @signature JstFinder(other);
 *
 * @param[in]   haystack    The haystack to be searched.
 * @param[in]   other       Other finder of the same type (copy constructor).
 */

/*!
 * @mfn JstFinder#RegisteredExtensionPoint
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Returns the type of the registered enxtension point.
 *
 * @signature RegisteredExtensionPoint<TFinder>::Type;
 * @tparam  TFinder The type of the finder to get the registered extension point for.
 *
 * @return TExtensionPoint  The type of the registered extension point. See @link FinderExtensionPoint @endlink.
 */

template <typename TContainer, typename TPattern, typename TSpec>
struct Finder_<TContainer, TPattern, Jst<TSpec> >
{
    typedef Finder_<TContainer, TPattern, Jst<TSpec> > TFinder;
    typedef typename RegisteredExtensionPoint<TFinder>::Type TFinderExtensionPoint;

    TContainer*             _containerPtr;
    TFinderExtensionPoint   _extensionFunctor;
    bool                    _needReinit;

    Finder_() : _containerPtr(NULL), _needReinit(true)
    {}


    Finder_(TContainer & hystk) : _containerPtr(&hystk), _extensionFunctor(), _needReinit(true)
    {}

    // Copy constructor.
    Finder_(Finder_ const & other) : _containerPtr(other._containerPtr),
                                     _extensionFunctor(other._extensionFunctor),
                                     _needReinit(other._needReinit)
    {}

    // Assignment Operator
    Finder_ & operator=(Finder_ const & other)
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
SEQAN_CONCEPT_IMPL((JstTraversalConcept), Finder_<TContainer, TPattern, Jst<TSpec> >);

template <typename TContainer, typename TPattern, typename TSpec>
SEQAN_CONCEPT_IMPL((JstTraversalConcept), Finder_<TContainer, TPattern, Jst<TSpec> > const);

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct ContextIteratorPosition<Finder_<TContainer, TPattern, TSpec> >
{
    typedef Finder_<TContainer, TPattern, TSpec> TFinder_;
    typedef typename RegisteredExtensionPoint<TFinder_>::Type TFinderFunctor_;
    typedef typename ContextIteratorPosition<TFinderFunctor_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct ContextIteratorPosition<Finder_<TContainer, TPattern, TSpec> const > :
    ContextIteratorPosition<Finder_<TContainer, TPattern, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct RequireFullContext<Finder_<TContainer, TPattern, TSpec> >
{
    typedef Finder_<TContainer, TPattern, TSpec> TFinder_;
    typedef typename RegisteredExtensionPoint<TFinder_>::Type TFinderFunctor_;
    typedef typename RequireFullContext<TFinderFunctor_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct RequireFullContext<Finder_<TContainer, TPattern, TSpec> const > :
    RequireFullContext<Finder_<TContainer, TPattern, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct Container<Finder_<TContainer, TPattern, Jst<TSpec> > >
{
    typedef TContainer Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct Container<Finder_<TContainer, TPattern, Jst<TSpec> > const>
{
    typedef TContainer const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct GetState<Finder_<TContainer, TPattern, Jst<TSpec> > >
{
    typedef Finder_<TContainer, TPattern, Jst<TSpec> > TFinder_;
    typedef typename RegisteredExtensionPoint<TFinder_>::Type TFinderExtension_;
    typedef typename GetState<TFinderExtension_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct GetState<Finder_<TContainer, TPattern, Jst<TSpec> > const> :
    GetState<Finder_<TContainer, TPattern, Jst<TSpec> > >{};

// ----------------------------------------------------------------------------
// Metafunction GetJstTraverser
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct GetJstTraverser<Finder_<TContainer, TPattern, TSpec> >
{
    typedef Finder_<TContainer, TPattern, TSpec> TFinder_;
    typedef typename ContextIteratorPosition<TFinder_>::Type TContextPosition_;
    typedef typename RequireFullContext<TFinder_>::Type TRequireContext_;
    typedef typename GetState<TFinder_>::Type TState;
    typedef JstTraverser<TContainer, TState, JstTraverserSpec<TContextPosition_, TRequireContext_> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
inline typename GetState<Finder_<TContainer, TPattern, Jst<TSpec> > const>::Type
getState(Finder_<TContainer, TPattern, Jst<TSpec> > const & finder)
{
    return getState(finder._extensionFunctor);  // Delegates to extension functor.
}

// ----------------------------------------------------------------------------
// Function setState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TState>
inline void
setState(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
         TState const & state)
{
    setState(finder._extensionFunctor, state);  // Delegates to extension functor.
}

// ----------------------------------------------------------------------------
// Function initState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
inline void
initState(Finder_<TContainer, TPattern, Jst<TSpec> > & finder)
{
    initState(finder._extensionFunctor);  // Delegates to extension functor.
}

// ----------------------------------------------------------------------------
// Function execute()
// ----------------------------------------------------------------------------

template <typename TResult, typename TFinderExtension, typename TContextIter>
inline void
execute(TResult & res,
        TFinderExtension & extensionFunctor,
        TContextIter & contextIter)
{
    extensionFunctor(res, contextIter);
}

// ----------------------------------------------------------------------------
// Function deliverContext()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate, typename TTraverser, typename TTag>
inline typename Size<TTraverser>::Type
deliverContext(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
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
// Function container()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
inline typename Container<Finder_<TContainer, TPattern, Jst<TSpec> > >::Type &
container(Finder_<TContainer, TPattern, Jst<TSpec> > & finder)
{
    return *finder._containerPtr;
}

template <typename TContainer, typename TPattern, typename TSpec>
inline typename Container<Finder_<TContainer, TPattern, Jst<TSpec> > const>::Type &
container(Finder_<TContainer, TPattern, Jst<TSpec> > const & finder)
{
    return *finder._containerPtr;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TExtensionSpec, typename TPattern, typename TScoreLimit>
inline void
init(FinderExtensionPoint<TFinder, TExtensionSpec> & extensionFunctor,
     TPattern & pattern,
     TScoreLimit const & /*scoreLimit*/)
{
    init(extensionFunctor, pattern);
}

//// ----------------------------------------------------------------------------
//// Function _find()
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
//inline void
//_find(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
//      TPattern & pattern,
//      TDelegate & delegate,
//      int scoreLimit,
//      Tag<T>)
//{
//}

// ----------------------------------------------------------------------------
// Function find()                                                     [Serial]
// ----------------------------------------------------------------------------

/*!
 * @fn JstFinder#find
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Triggers the search over the haystack.
 *
 * @signature find(finder, pattern, delegate[, limit]);
 *
 * @param[in,out]   finder      The finder which manages the search.
 * @param[in,out]   pattern     The pattern to be searched. Of type @link Pattern @endlink.
 * @param[in,out]   delegate    An additional functor called on success. See @link JstFinderExtensionConcept#execute @endlink.
 * @param[in]       limit       An optional parameter setting the score limit (<tt><= 0</tt>).
 */

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate, typename TParallel>
inline void
find(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
     TPattern & pattern,
     TDelegate & delegate,
     int scoreLimit,
     Tag<TParallel> const& parallelTag)
{
    typedef Finder_<TContainer, TPattern, Jst<TSpec> > TFinder;
    typedef typename GetJstTraverser<TFinder>::Type TTraverser;

    // Set up the journaled string tree traversal.
    TTraverser traverser(container(finder), length(needle(pattern)) - scoreLimit);

    if (finder._needReinit)
    {
        init(finder._extensionFunctor, pattern, scoreLimit);
        finder._needReinit = false;
    }
    traverse(finder, delegate, traverser, parallelTag);
}

template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
inline void
find(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
     TPattern & pattern,
     TDelegate & delegate,
     int scoreLimit = 0)
{
    find(finder, pattern, delegate, scoreLimit, Serial());
}

//template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
//inline void
//find(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
//     TPattern & pattern,
//     TDelegate & delegate,
//     int scoreLimit)
//{
//    find(finder, pattern, delegate, scoreLimit, 1);
//}
//
//template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
//inline void
//find(Finder_<TContainer, TPattern, Jst<TSpec> > & finder,
//     TPattern & pattern,
//     TDelegate & delegate)
//{
//    find(finder, pattern, delegate, 0);
//}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_IMPL_H_
