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

// TODO(rmaerker): Finder2 class needs to be documented.
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
 *            struct Finder2<THaystack, TPattern, JstFinder>;
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
struct Finder2<TContainer, TPattern, Jst<TSpec> >
{
    typedef Finder2<TContainer, TPattern, Jst<TSpec> > TFinder;
    typedef typename RegisteredExtensionPoint<TFinder>::Type TFinderExtensionPoint;

    TContainer*             _containerPtr;
    TFinderExtensionPoint   _extensionFunctor;
    Pair<unsigned>          _traversalParams;

    Finder2() : _containerPtr(NULL), _traversalParams(Pair<unsigned>(0, 0))
    {}


    Finder2(TContainer & hystk) : _containerPtr(&hystk), _extensionFunctor(), _traversalParams(Pair<unsigned>(0, 0))
    {}

    // Copy constructor.
    Finder2(Finder2 const & other) : _containerPtr(other._containerPtr),
                                     _extensionFunctor(other._extensionFunctor),
                                     _traversalParams(other._traversalParams)

    {}

    // Assignment Operator
    Finder2 & operator=(Finder2 const & other)
    {
        if (this != &other)
        {
            _containerPtr = other._containerPtr;
            _extensionFunctor = other._extensionFunctor;
            _traversalParams = other._traversalParams;
        }
        return *this;
    }
};

template <typename TContainer, typename TPattern, typename TSpec>
SEQAN_CONCEPT_IMPL((JstTraversalConcept), Finder2<TContainer, TPattern, Jst<TSpec> >);

template <typename TContainer, typename TPattern, typename TSpec>
SEQAN_CONCEPT_IMPL((JstTraversalConcept), Finder2<TContainer, TPattern, Jst<TSpec> > const);

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
    typedef typename RegisteredExtensionPoint<TFinder_>::Type TFinderFunctor_;
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
    typedef typename RegisteredExtensionPoint<TFinder_>::Type TFinderFunctor_;
    typedef typename RequireFullContext<TFinderFunctor_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct RequireFullContext<Finder2<TContainer, TPattern, TSpec> const > :
    RequireFullContext<Finder2<TContainer, TPattern, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct Container<Finder2<TContainer, TPattern, Jst<TSpec> > >
{
    typedef TContainer Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct Container<Finder2<TContainer, TPattern, Jst<TSpec> > const>
{
    typedef TContainer const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
struct GetState<Finder2<TContainer, TPattern, Jst<TSpec> > >
{
    typedef Finder2<TContainer, TPattern, Jst<TSpec> > TFinder_;
    typedef typename RegisteredExtensionPoint<TFinder_>::Type TFinderExtension_;
    typedef typename GetState<TFinderExtension_>::Type Type;
};

template <typename TContainer, typename TPattern, typename TSpec>
struct GetState<Finder2<TContainer, TPattern, Jst<TSpec> > const> :
    GetState<Finder2<TContainer, TPattern, Jst<TSpec> > >{};

// ----------------------------------------------------------------------------
// Metafunction GetJstTraverser
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
// Function getGetState()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec>
inline typename GetState<FinderExtensionPoint<TFinder, TSpec> >::Type
getState(FinderExtensionPoint<TFinder, TSpec> const & /*extensionFunctor*/)
{
    return typename GetState<FinderExtensionPoint<TFinder, TSpec> >::Type();
}

// ----------------------------------------------------------------------------
// Function setGetState()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TSpec, typename TState>
inline void
setState(FinderExtensionPoint<TFinder, TSpec> const & /*extensionFunctor*/,
         TState const & /*state*/)
{
    // no-op function.
}

// ----------------------------------------------------------------------------
// Function getState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec>
inline typename GetState<Finder2<TContainer, TPattern, Jst<TSpec> > const>::Type
getState(Finder2<TContainer, TPattern, Jst<TSpec> > const & finder)
{
    return getState(finder._extensionFunctor);  // Delegates to extension functor.
}

// ----------------------------------------------------------------------------
// Function setState()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TState>
inline void
setState(Finder2<TContainer, TPattern, Jst<TSpec> > & finder,
         TState const & state)
{
    setState(finder._extensionFunctor, state);  // Delegates to extension functor.
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
deliverContext(Finder2<TContainer, TPattern, Jst<TSpec> > & finder,
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
inline typename Container<Finder2<TContainer, TPattern, Jst<TSpec> > >::Type &
container(Finder2<TContainer, TPattern, Jst<TSpec> > & finder)
{
    return *finder._containerPtr;
}

template <typename TContainer, typename TPattern, typename TSpec>
inline typename Container<Finder2<TContainer, TPattern, Jst<TSpec> > const>::Type &
container(Finder2<TContainer, TPattern, Jst<TSpec> > const & finder)
{
    return *finder._containerPtr;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TFinder, typename TExtensionSpec, typename TPattern, typename TScoreLimit>
inline Pair<unsigned>
init(FinderExtensionPoint<TFinder, TExtensionSpec> & extensionFunctor,
     TPattern & pattern,
     TScoreLimit const & /*scoreLimit*/)
{
    return init(extensionFunctor, pattern);
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TPattern, typename TSpec, typename TState, typename TErrors>
inline void
init(Finder2<TContainer, TPattern, Jst<TSpec> > & finder,
     TPattern & pattern,
     FinderState<TState> & state,
     TErrors const & errors)
{
    finder._traversalParams = init(finder.extensionFunctor, pattern, state, errors);
}

template <typename TContainer, typename TPattern, typename TSpec, typename TErrors>
inline void
init(Finder2<TContainer, TPattern, Jst<TSpec> > & finder,
     TPattern & pattern,
     FinderState<Nothing> & /*state*/,
     TErrors const & errors)
{
    finder._traversalParams = init(finder.extensionFunctor, pattern, errors);
}

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

//template <typename TContainer, typename TPattern, typename TSpec, typename TDelegate>
//inline void
//find(Finder2<TContainer, TPattern, Jst<TSpec> > & finder,
//     TPattern & pattern,
//     TDelegate & delegate,
//     int scoreLimit = 0)
//{
//    typedef Finder2<TContainer, TPattern, Jst<TSpec> > TFinder;
//    typedef typename GetJstTraverser<TFinder>::Type TTraverser;
//
//    // Set up the journaled string tree traversal.
//    TTraverser traverser(container(finder), length(needle(pattern)) - scoreLimit);
//#ifdef PROFILE_JST_INTERN
//    double timeBuild = sysTime();
//#endif
//    while (journalNextBlock(container(finder), contextSize(traverser)))  // Generating the journal index for the next block.
//    {
//        if (finder._needReinit)
//        {
//            init(finder._extensionFunctor, pattern, scoreLimit);
//            finder._needReinit = false;
//        }
//#ifdef PROFILE_JST_INTERN
//        std::cerr << "Time-C: " << sysTime() - timeBuild << " s" << std::endl;
//        double timeSearch = sysTime();
//#endif
//        traverse(finder, delegate, traverser);
//#ifdef PROFILE_JST_INTERN
//        std::cerr << "Time-S: " << sysTime() - timeSearch << " s" << std::endl;
//        timeBuild = sysTime();
//#endif
//    }
//}

template <typename THystk, typename TPattern, typename TDelegate, typename TState, typename TErrors, typename TSpec>
void find(THystk & jst,
          TPattern & pattern,
          TDelegate & delegate,
          FinderState<TState> & state,
          TErrors const & errors,
          Jst<TSpec> const & /*tag*/)
{
    typedef Finder2<THystk, TPattern, Jst<TSpec> > TFinder;
    typedef typename GetJstTraverser<TFinder>::Type TTraverser;

    TFinder finder(jst);
    init(finder, pattern, state, errors);
    TTraverser traverser(container(finder), finder._traversalParams.i1, finder._traversalParams.i2);
    traverse(finder, delegate, traverser);
}

// Without errors
template <typename THystk, typename TPattern, typename TDelegate, typename TState, typename TSpec>
void find(THystk & jst,
          TPattern & pattern,
          TDelegate & delegate,
          FinderState<TState> & state,
          Jst<TSpec> const & tag)
{
    find(jst, pattern, delegate, state, 0, tag);
}

// Without state.
template <typename THystk, typename TPattern, typename TDelegate, typename TErrors, typename TSpec>
void find(THystk & jst,
          TPattern & pattern,
          TDelegate & delegate,
          TErrors const & errors,
          Jst<TSpec> const & tag)
{
    FinderState<> defaultState;
    find(jst, pattern, delegate, defaultState, errors, tag);
}

// Without state and errors.
template <typename THystk, typename TPattern, typename TDelegate, typename TSpec>
void find(THystk & jst,
          TPattern & pattern,
          TDelegate & delegate,
          Jst<TSpec> const & tag)
{
    find(jst, pattern, delegate, 0, tag);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FIND_JOURNALED_STRING_TREE_IMPL_H_
