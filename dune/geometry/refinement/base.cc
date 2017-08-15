// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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

    //! Get the number of Vertices
    static int nVertices(int level, Dune::VirtualRefinementTag::Level = {});
    static int nVertices(int nIntervals, Dune::VirtualRefinementTag::PerAxis);
    //! Get a VertexIterator
    static VertexIterator vBegin(int level, Dune::VirtualRefinementTag::Level = {});
    static VertexIterator vBegin(int nIntervals, Dune::VirtualRefinementTag::PerAxis);
    //! Get a VertexIterator
    static VertexIterator vEnd(int level, Dune::VirtualRefinementTag::Level = {});
    static VertexIterator vEnd(int nIntervals, Dune::VirtualRefinementTag::PerAxis);

    //! Get the number of Elements
    static int nElements(int level, Dune::VirtualRefinementTag::Level = {});
    static int nElements(int nIntervals, Dune::VirtualRefinementTag::PerAxis);
    //! Get an ElementIterator
    static ElementIterator eBegin(int level, Dune::VirtualRefinementTag::Level = {});
    static ElementIterator eBegin(int nIntervals, Dune::VirtualRefinementTag::PerAxis);
    //! Get an ElementIterator
    static ElementIterator eEnd(int level, Dune::VirtualRefinementTag::Level = {});
    static ElementIterator eEnd(int nIntervals, Dune::VirtualRefinementTag::PerAxis);
  };

  /*! \} */

  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  int StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::nVertices(int level, Dune::VirtualRefinementTag::Level)
  {
    return RefinementImp::nVertices(1<<level);
  }
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  int StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::nVertices(int nIntervals, Dune::VirtualRefinementTag::PerAxis)
  {
    return RefinementImp::nVertices(nIntervals);
  }

  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::VertexIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::vBegin(int level, Dune::VirtualRefinementTag::Level)
  {
    return RefinementImp::vBegin(1<<level);
  }
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::VertexIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::vBegin(int nIntervals, Dune::VirtualRefinementTag::PerAxis)
  {
    return RefinementImp::vBegin(nIntervals);
  }

  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::VertexIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::vEnd(int level, Dune::VirtualRefinementTag::Level)
  {
    return RefinementImp::vEnd(1<<level);
  }
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::VertexIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::vEnd(int nIntervals, Dune::VirtualRefinementTag::PerAxis)
  {
    return RefinementImp::vEnd(nIntervals);
  }

  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  int StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::nElements(int level, Dune::VirtualRefinementTag::Level)
  {
    return RefinementImp::nElements(1<<level);
  }
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  int StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::nElements(int nIntervals, Dune::VirtualRefinementTag::PerAxis)
  {
    return RefinementImp::nElements(nIntervals);
  }

  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::ElementIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::eBegin(int level, Dune::VirtualRefinementTag::Level)
  {
    return RefinementImp::eBegin(1<<level);
  }
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::ElementIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::eBegin(int nIntervals, Dune::VirtualRefinementTag::PerAxis)
  {
    return RefinementImp::eBegin(nIntervals);
  }

  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::ElementIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::eEnd(int level, Dune::VirtualRefinementTag::Level)
  {
    return RefinementImp::eEnd(1<<level);
  }
  template<unsigned topologyId, class CoordType,
      unsigned coerceToId, int dimension_>
  typename StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::ElementIterator
  StaticRefinement<topologyId, CoordType, coerceToId, dimension_>::eEnd(int nIntervals, Dune::VirtualRefinementTag::PerAxis)
  {
    return RefinementImp::eEnd(nIntervals);
  }


} // namespace Dune

#endif // DUNE_GEOMETRY_REFINEMENT_BASE_CC
