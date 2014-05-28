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
// The concept for the finder extension.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FINDER_EXTENSION_POINT_CONCEPT_H_
#define EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FINDER_EXTENSION_POINT_CONCEPT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename TExtensionPoint>
struct ExtensionRegistry_;

template <typename T>
struct GetPattern;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @concept JstFinderExtensionConcept
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Defines the interface to extend the @link JstFinder @endlink by a customized algorithm.
 *
 * @signature concept JstFinderExtensionConcept;
 *
 * This concept provides an interface to generically extend the @link JstFinder @endlink with customized
 * algorithms. The new algorithm must be implemented in form of an functor of type @link FinderExtensionPoint @endlink.
 *
 */

/*!
 * @mfn JstFinderExtensionConcept#GetState
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Returns the type of the state used for the external algorithm.
 *
 *  Note this metafunction is optional and defaults to @link Nothing @endlink if not overloaded.
 *
 * @signature GetState<TExtension>::Type;
 * @tparam TExtension Type of the @link FinderExtensionPoint @endlink.
 *
 * @return TState The type of the state used by <tt>TExtension</tt>.
 */

/*!
 * @mfn JstFinderExtensionConcept#RegisteredExtensionPoint
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Returns the type of the external algorithm registered to the @link JstFinder @endlink.
 *
 * @signature RegisteredExtensionPoint<TFinder>::Type;
 * @tparam TFinder The type of the finder for which the type of the extension point should be determined.
 *
 * @return TExtension The type of the registered extension point.
 */

/*!
 * @mfn JstFinderExtensionConcept#ContextIteratorPosition
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Returns the specialization for the iterator position within the context.
 *
 * Note this metafunction is optional and defaults to @link JstTraverserContextPositionTags#ContextPositionLeft @endlink
 * if not overloaded. See @link JstTraverser @endlink for more information.
 *
 * @signature ContextIteratorPosition<TExtension>::Type
 * @tparam TExtension The type of the extension point for which the context position should be specialized.
 *         Must only be one of @link JstTraverserContextPositionTags @endlink
 *
 *  @return TContextPosition The tag specializing the context positio of the iterator.
 */

/*!
 * @mfn JstFinderExtensionConcept#RequireFullContext
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Returns a logical type to indicate whether the algorithm requires the full context to be accessible.
 *
 * Note this metafunction is optional and defaults to @link LogicalValuesTags#True @endlink
 * See @link JstTraverser @endlink for more information.
 *
 * @signature RequireFullContext<TExtension>::Type
 * @tparam TExtension The type of the extension point for which the context should be specialized.
 *         Must only be one of @link LogicalValuesTags @endlink
 *
 *  @return TLogic @link LogicalValuesTags#True @endlink if the entire context needs to be accessible,
 *  @link LogicalValuesTags#False @endlink otherwise.
 */

/*!
 * @fn JstFinderExtensionConcept#execute
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Executes the algorithm registered to the @link JstFinder @endlink.
 *
 * Note this function is optional and executes the default functor @link FinderExtensionPoint @endlink. Overload this
 * function if a special behavior is needed to execute the extension.
 *
 * @signature execute(res, ext, it);
 *
 * @param[in,out] res   The result value of type @link Pair @endlink with the first value being <tt>bool</tt> indicating
 * a success event if <tt>true</tt> and the second value being of type <tt>size_t</tt> to return the step size. Per
 * default <tt>res.i1 = false </tt> and <tt>res.i2 = 1 </tt>.
 *
 * @param[in,out] ext   The functor of type @link FinderExtensionPoint @endlink to be executed.
 * @param[in]     it    The iterator to the current context.
 */

/*!
 * @fn JstFinderExtensionConcept#init
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Initializes the external algorithm registered to the @link JstFinder @endlink.
 *
 * @signature init(ext, obj);
 * @signature init(ext, obj, limit);
 *
 * @param[in,out]   ext     The functor of type @link FinderExtensionPoint @endlink to be initialized.
 * @param[in]       obj     The object usued to initialize the @link FinderExtensionPoint @endlink with.
 * @param[in]       limit   The score limit. Must only be less or equal than <tt>0</tt>.
 */

/*!
 * @fn JstFinderExtensionConcept#getState
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Returns the state of the external algorithm registered to the @link JstFinder @endlink.
 *
 * Note this function is optional and returns @link Nothing @endlink if no state is needed. Overload this function
 * to return the state for a specialized @link FinderExtensionPoint @endlink.
 *
 * @signature TState getState(ext);
 *
 * @param[in]   ext The functor of type @link FinderExtensionPoint @endlink to get the state for.
 *
 * @return TState The state of the @link FinderExtensionPoint @endlink of type
 *                @link JstFinderExtensionConcept#GetState @endlink.
 *
 * @see JstFinderExtensionConcept#setState
 */

/*!
 * @fn JstFinderExtensionConcept#setState
 * @headerfile <seqan/find_journaled_string_tree.h>
 * @brief Sets the state of the external algorithm registered to the @link JstFinder @endlink.
 *
 * Note this function is optional and defaults to a non-operational function. Overload this function
 * to set the state for a specialized @link FinderExtensionPoint @endlink.
 *
 * @signature TState setState(ext, state);
 *
 * @param[in, out]   ext    The functor of type @link FinderExtensionPoint @endlink to set the state for.
 * @param[in]        state  The new state to be set.
 *
 * @see JstFinderExtensionConcept#getState
 */

template <typename TExtension>
struct JstFinderExtensionConcept :
    DefaultConstructible<TExtension>
{
    typedef typename GetState<TExtension>::Type TState;
    typedef typename ExtensionRegistry_<TExtension>::Type TExtensionRegistry;
    typedef typename GetPattern<TExtensionRegistry>::Type TPattern;
    typedef typename Iterator<CharString, Standard>::Type TIter;
    typedef typename Container<TExtensionRegistry>::Type  TContainer;
    typedef typename Size<TContainer>::Type               TSize;

    SEQAN_CONCEPT_ASSERT((JstTraversalConcept<TExtensionRegistry>));

    TExtension        e;
    TState            state;
    TPattern          pattern;
    int               limit;
    TIter             it;
    Pair<bool, TSize> res;


    SEQAN_CONCEPT_USAGE(JstFinderExtensionConcept)
    {
        sameType(getState(e), state);
        state = getState(e);
        setState(e, state);

        execute(res, e, it);

        init(e, pattern, limit);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}

#endif // EXTRAS_INCLUDE_SEQAN_FIND_JOURNALED_STRING_TREE_FINDER_EXTENSION_POINT_CONCEPT_H_
