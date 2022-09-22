// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_REFINEMENT_BASE_CC
#define DUNE_GEOMETRY_REFINEMENT_BASE_CC

/*!
 * \file
 *
 * \brief This file contains the parts independent of a particular
 *        \ref Refinement implementation.
 */

#include <dune/geometry/type.hh>

namespace Dune
{
  /*!
   * \addtogroup Refinement Refinement
   * \{
   */

  /*!
   * \brief This namespace contains the implementation of \ref
   *        Refinement.
   */
  namespace RefinementImp
  {
    // /////////////////////////////////
    //
    // Declaration of RefinementImp::Traits
    //

#ifdef DOXYGEN
    // This is just for Doxygen
    /*!
     * \brief Mapping from \a geometryType, \a CoordType and \a coerceTo
     *        to a particular \ref Refinement implementation.
     *
     * \tparam topologyId The topology id of the element to refine
     * \tparam CoordType  The C++ type of the coordinates
     * \tparam coerceToId The topologyId of the subelements
     * \tparam dimension  The dimension of the refinement.
     * \tparam Dummy      Dummy parameter which can be used for SFINAE, should
     *                    always be void.
     *
     * Each \ref Refinement implementation has to define one or more
     * specialisations of this struct to declare what it implements.
     * Template class Refinement uses this struct to know which
     * implementation it should inherit from.  Since non-type template
     * arguments of specializations may not involve template parameters, it is
     * often impossible to specify the specialization for all cases directly.
     * As the workaround, the template parameter \a Dummy can be used for
     * SFINAE with \a enable_if.
     *
     * Each specialisation should contain a single member typedef Imp,
     * e.g.:
     * \code
     * template<class CoordType>
     * struct Traits<sphereTopologyId, CoordType, Impl::CubeToplogy<2>::id, 2>
     * {
     *   typedef SquaringTheCircle::Refinement Imp;
     * };
     * \endcode
     */
    template<unsigned topologyId, class CoordType,
        unsigned coerceToId, int dimension, class Dummy = void>
    struct Traits
    {
      //! The implementation this specialisation maps to
      typedef SquaringTheCircle::Refinement Imp;
    };

#else // !DOXYGEN

    // Doxygen won't see this

    template<unsigned topologyId, class CoordType,
        unsigned coerceToId, int dimension, class = void>
    struct Traits;

#endif // !DOXYGEN
  } // namespace RefinementImp


  /*!
   * \brief Holds the number of refined intervals per axis needed for virtual and static refinement.
   *
   * To create an object of this class, call either refinementIntervals() or refinementLevels(). The first on will just
   * pass its input to the constructor, the latter one will pass 2^{input} to the constructor to be consistent with the
   * meaning of levels in a grid context.
   */
  class RefinementIntervals{
    int intervals_=1;

  public:
    explicit RefinementIntervals(int i) : intervals_(i) {}

    int intervals() const { return intervals_; }
  };

  /*!
   * \brief Creates a RefinementIntervals object
   *
   * \param intervals Number of refined intervals per axis
   */
  inline RefinementIntervals refinementIntervals(int intervals)
  {
    return RefinementIntervals{intervals};
  }
  /*!
   * \brief Creates a RefinementIntervals object
   *
   * \param levels Number of refinement levels, translates to \f$2^{levels}\f$ intervals per axis
   */
  inline RefinementIntervals refinementLevels(int levels)
  {
    return RefinementIntervals{1<<levels};
  }

  // ///////////////
  //
  //  Static Refinement
  //

  /*!
   * \brief Wrap each \ref Refinement implementation to get a
   *        consistent interface
   *
   * \tparam topologyId The topology id of the element to refine
   * \tparam CoordType  The C++ type of the coordinates
   * \tparam coerceToId The topology id of the subelements
   * \tparam dimension  The dimension of the refinement.
   */
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  class StaticRefinement
    : public RefinementImp::Traits<topologyId, CoordType,
          coerceToId, dimension_ >::Imp
  {
  public:
#ifdef DOXYGEN
    /*!
     * \brief The Codim struct inherited from the \ref Refinement implementation
     *
     * \tparam codimension There is a different struct Codim for each codimension
     */
    template<int codimension>
    struct Codim
    {
      /*!
       * \brief The SubEntityIterator for each codim
       *
       * This is \em some sort of type, not necessarily a typedef
       */
      typedef SubEntityIterator;
    };

    //! The VertexIterator of the Refinement
    typedef Codim<dimension>::SubEntityIterator VertexIterator;
    //! The ElementIterator of the Refinement
    typedef Codim<0>::SubEntityIterator ElementIterator;

    /*!
     * \brief The CoordVector of the Refinement
     *
     * This is always a typedef to a FieldVector
     */
    typedef CoordVector;

    /*!
     * \brief The IndexVector of the Refinement
     *
     * This is always a typedef to a FieldVector
     */
    typedef IndexVector;
#endif

    typedef typename RefinementImp::Traits< topologyId, CoordType, coerceToId, dimension_>::Imp RefinementImp;

    using RefinementImp::dimension;

    using RefinementImp::Codim;

    using typename RefinementImp::VertexIterator;
    using typename RefinementImp::CoordVector;

    using typename RefinementImp::ElementIterator;
    using typename RefinementImp::IndexVector;

    /*!
     * \brief Get the number of Vertices
     *
     * \param tag RefinementIntervals object returned by either refinementIntervals() or refinementLevels()
     */
    static int nVertices(Dune::RefinementIntervals tag)
    {
      return RefinementImp::nVertices(tag.intervals());
    }

    /*!
     * \brief Get a VertexIterator
     *
     * \param tag RefinementIntervals object returned by either refinementIntervals() or refinementLevels()
     */
    static VertexIterator vBegin(Dune::RefinementIntervals tag)
    {
      return RefinementImp::vBegin(tag.intervals());
    }

    /*!
     * \brief Get a VertexIterator
     *
     * \param tag RefinementIntervals object returned by either refinementIntervals() or refinementLevels()
     */
    static VertexIterator vEnd(Dune::RefinementIntervals tag)
    {
      return RefinementImp::vEnd(tag.intervals());
    }

    /*!
     * \brief Get the number of Elements
     *
     * \param tag RefinementIntervals object returned by either refinementIntervals() or refinementLevels()
     */
    static int nElements(Dune::RefinementIntervals tag)
    {
      return RefinementImp::nElements(tag.intervals());
    }

    /*!
     * \brief Get an ElementIterator
     *
     * \param tag RefinementIntervals object returned by either refinementIntervals() or refinementLevels()
     */
    static ElementIterator eBegin(Dune::RefinementIntervals tag)
    {
      return RefinementImp::eBegin(tag.intervals());
    }

    /*!
     * \brief Get an ElementIterator
     *
     * \param tag RefinementIntervals object returned by either refinementIntervals() or refinementLevels()
     */
    static ElementIterator eEnd(Dune::RefinementIntervals tag)
    {
      return RefinementImp::eEnd(tag.intervals());
    }
  };

  /*! \} */
} // namespace Dune

#endif // DUNE_GEOMETRY_REFINEMENT_BASE_CC
