// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_REFINEMENT_HCUBE_CC
#define DUNE_GEOMETRY_REFINEMENT_HCUBE_CC

/*!
 * \file
 * \brief This file contains the \ref Refinement implementation for
 *        hypercubes (quadrilaterals, hexahedrons, etc.).
 *
 * See \ref HCubeRefinement.
 */

/*!
 * \defgroup HCubeRefinement Refinement implementation for hypercubes
 *  \ingroup Refinement
 *
 * \section Iterators The Iterators
 * <!--=========================-->
 *
 * For the iterators we have to hack around a bit.  The problem is as
 * follows:
 * \code
 * template<int A>
 * class outer
 * {
 *   template<int B>
 *   class inner;
 * };
 * \endcode
 * C++ does not allow specialisation of the inner class when the outer
 * class is not specialized.
 *
 * So I had to create a baseclass for the iterators which is not inside
 * another class.  This base class can then be specialized, and the
 * real Iterator class inherits from it.  I gave it the somewhat clumsy
 * name RefinementSubEntityIteratorSpecial.
 */

#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/iteratorfacades.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

#include "base.cc" // for RefinementTraits

namespace Dune
{
  namespace RefinementImp
  {
    /*!
     * \brief This namespace contains the \ref Refinement implementation
     * for hypercubes (GeometryType::cube).
     *
     * See \ref HCubeRefinement.
     */
    namespace HCube
    {
      /*!
       * \brief \ref Refinement implementation for hypercubes
       *
       * \param dimension_ Dimension of the refined hypercube
       * \param CoordType  Coordinate type of the refined hypercube
       *
       *  The interface is the same as for \ref Dune::StaticRefinement (apart
       * from the template parameters).
       */
      template<int dimension_, class CoordType>
      class RefinementImp
      {
      public:
        /** \brief Know your own dimension \hideinitializer */
        constexpr static int dimension = dimension_;
        //- Know yourself
        typedef RefinementImp<dimension, CoordType> Refinement;

        template<int codimension>
        struct Codim;
        typedef typename Codim<dimension>::SubEntityIterator VertexIterator;
        typedef FieldVector<CoordType, dimension> CoordVector;
        typedef typename Codim<0>::SubEntityIterator ElementIterator;
        typedef FieldVector<int, (1<<dimension)> IndexVector;

        static unsigned nVertices(unsigned nIntervals);
        static VertexIterator vBegin(unsigned nIntervals);
        static VertexIterator vEnd(unsigned nIntervals);

        static unsigned nElements(unsigned nIntervals);
        static ElementIterator eBegin(unsigned nIntervals);
        static ElementIterator eEnd(unsigned nIntervals);
      };

      template<int dimension, class CoordType>
      template<int codimension>
      struct RefinementImp<dimension, CoordType>::Codim
      {
        class SubEntityIterator;
        typedef Dune::AxisAlignedCubeGeometry<CoordType,dimension-codimension,dimension> Geometry;
      };

      template<int dimension, class CoordType>
      unsigned
      RefinementImp<dimension, CoordType>::
      nVertices(unsigned nIntervals)
      {
        // return (nIntervals + 1)^dim
        return Dune::power(nIntervals+1u, unsigned(dimension));
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vBegin(unsigned nIntervals)
      {
        return VertexIterator(0,nIntervals);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vEnd(unsigned nIntervals)
      {
        return VertexIterator(nVertices(nIntervals),nIntervals);
      }

      template<int dimension, class CoordType>
      unsigned
      RefinementImp<dimension, CoordType>::
      nElements(unsigned nIntervals)
      {
        static_assert(dimension >= 0,
                      "Negative dimension given, what the heck is that supposed to mean?");
        // return nIntervals^dim
        return Dune::power(nIntervals, unsigned(dimension));
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eBegin(unsigned nIntervals)
      {
        return ElementIterator(0,nIntervals);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eEnd(unsigned nIntervals)
      {
        return ElementIterator(nElements(nIntervals),nIntervals);
      }

      //
      // The iterators
      //

#ifdef DOXYGEN
      /*!
       * \brief SubEntityIterator base class for hypercube refinement
       *
       * \tparam dimension   Dimension of the refined element
       * \tparam CoordType   Coordinate type of the refined element
       * \tparam codimension Codimension of the iterator
       *
       * This is the base class for SubEntityIterators.  We have to use
       * this construct because RefinementImp<...>::%codim<...> cannot
       * be specialized without first specializing RefinementImp.
       */
      template<int dimension, class CoordType, int codimension>
      class RefinementSubEntityIteratorSpecial {};
#else //!DOXYGEN
      template<int dimension, class CoordType, int codimension>
      class RefinementSubEntityIteratorSpecial;
#endif //DOXYGEN

      // for vertices

      template<int dimension, class CoordType>
      class RefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::template Codim<dimension>::SubEntityIterator Common;
        typedef typename Refinement::CoordVector CoordVector;

        CoordVector coords() const;

      private:
        const Common & asCommon() const
        {
          return *static_cast<const Common*>(this);
        }
      };

      template<int dimension, class CoordType>
      typename RefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>::CoordVector
      RefinementSubEntityIteratorSpecial<dimension, CoordType, dimension>::
      coords() const
      {
        std::array<unsigned int, dimension> v(asCommon().vertexCoord());
        CoordVector c;
        for (int d = 0; d < dimension; d++)
        {
          c[d] = v[d]*1.0 / asCommon()._nIntervals;
        }
        return c;
      }

      // for elements

      template<int dimension, class CoordType>
      class RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::template Codim<0>::SubEntityIterator Common;
        typedef typename Refinement::IndexVector IndexVector;
        typedef typename Refinement::CoordVector CoordVector;

        IndexVector vertexIndices() const;
        CoordVector coords() const;

      private:
        const Common & asCommon() const
        {
          return *static_cast<const Common*>(this);
        }
      };

      template<int dimension, class CoordType>
      typename RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::IndexVector
      RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::
      vertexIndices() const
      {
        constexpr static int nIndices = 1 << dimension;

        // cell index tuple
        std::array<unsigned int, dimension> e(asCommon().cellCoord());

        // vertices
        IndexVector vec;
        for(int i = 0; i < nIndices; ++i)
        {
          int base = 1;
          std::array<unsigned int, dimension> alpha(asCommon().idx2multiidx(i));
          for (int d = 0; d < dimension; d++) {
            vec[i] += (alpha[d] + e[d]) * base;
            base *= asCommon()._nIntervals+1;
          }
        }
        return vec;
      }

      template<int dimension, class CoordType>
      typename RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::CoordVector
      RefinementSubEntityIteratorSpecial<dimension, CoordType, 0>::
      coords() const
      {
        std::array<unsigned int, dimension> v(asCommon().cellCoord());
        CoordVector c;
        for (int d=0; d<dimension; d++)
        {
          c[d] = (v[d]*1.0 + 0.5) / asCommon()._nIntervals;
        }
        return c;
      }

      // common
      template<int dimension, class CoordType>
      template<int codimension>
      class RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator
        : public ForwardIteratorFacade<typename RefinementImp<dimension,
                  CoordType>::template Codim<codimension>::SubEntityIterator, int>,
          public RefinementSubEntityIteratorSpecial<dimension, CoordType, codimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::template Codim<codimension>::SubEntityIterator This;

        SubEntityIterator(unsigned int index, unsigned int nIntervals);

        bool equals(const This &other) const;
        void increment();

        int index() const;
        Geometry geometry () const;
      private:
        friend class RefinementSubEntityIteratorSpecial<dimension, CoordType, codimension>;
        unsigned int _index;
        unsigned int _nIntervals;

        std::array<unsigned int, dimension>
        cellCoord(unsigned int idx) const
        {
          return idx2coord(idx, _nIntervals);
        }

        std::array<unsigned int, dimension>
        vertexCoord(unsigned int idx) const
        {
          return idx2coord(idx, _nIntervals+1u);
        }

        std::array<unsigned int, dimension>
        cellCoord() const
        {
          return cellCoord(_index);
        }

        std::array<unsigned int, dimension>
        vertexCoord() const
        {
          return vertexCoord(_index);
        }

        std::array<unsigned int, dimension>
        idx2coord(unsigned int idx, unsigned int w) const
        {
          std::array<unsigned int, dimension> c;
          for (unsigned int d = 0; d < dimension; d++)
          {
            c[d] = idx%w;
            idx = idx/w;
          }
          return c;
        }

        unsigned int
        coord2idx(std::array<unsigned int, dimension> c, unsigned int w) const
        {
          unsigned int i = 0;
          for (unsigned int d = dimension; d > 0; d--)
          {
            i *= w;
            i += c[d-1];
          }
          return i;
        }

        std::array<unsigned int, dimension>
        idx2multiidx(unsigned int idx) const
        {
          std::array<unsigned int, dimension> alpha;
          for (unsigned int i = 0; i < dimension; ++i)
            alpha[i] = (idx >> i) & 1u;
          return alpha;
        }
      };

#ifndef DOXYGEN
      template<int dimension, class CoordType>
      template<int codimension>
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      SubEntityIterator(unsigned int index, unsigned int nIntervals)
        : _index(index), _nIntervals(nIntervals)
      {}

      template<int dimension, class CoordType>
      template<int codimension>
      bool
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      equals(const This &other) const
      {
        return ((_index == other._index) && (_nIntervals == other._nIntervals));
      }

      template<int dimension, class CoordType>
      template<int codimension>
      void
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      increment()
      {
        ++_index;
      }

      template<int dimension, class CoordType>
      template<int codimension>
      int
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      index() const
      {
        return _index;
      }

      template<int dimension, class CoordType>
      template<int codimension>
      typename RefinementImp<dimension, CoordType>::template Codim<codimension>::Geometry
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::geometry () const
      {
        std::array<unsigned int,dimension> intCoords = idx2coord(_index,_nIntervals);

        Dune::FieldVector<CoordType,dimension> lower;
        Dune::FieldVector<CoordType,dimension> upper;

        assert(codimension == 0 or codimension == dimension);

        if constexpr (codimension == 0) {
          for (size_t j = 0; j < dimension; j++)
          {
            lower[j] = double(intCoords[j])     / double(_nIntervals);
            upper[j] = double(intCoords[j] + 1) / double(_nIntervals);
          }

          return typename RefinementImp<dimension,
                                        CoordType>::template Codim<codimension>::Geometry(lower,upper);
        } else {
          for (size_t j = 0; j < dimension; j++)
            lower[j] = upper[j] = double(intCoords[j]) / double(_nIntervals);

          return typename RefinementImp<dimension,
                                        CoordType>::template Codim<codimension>::Geometry(lower,upper,std::bitset<dimension>(0));
        }
      }

#endif // DOXYGEN

    } // namespace HCube

    // ///////////////////////
    //
    // The refinement traits
    //

#ifndef DOXYGEN
    template<unsigned topologyId, class CoordType, unsigned coerceToId,
        int dim>
    struct Traits<
        topologyId, CoordType, coerceToId, dim,
        typename std::enable_if<
            (dim >= 2 &&
             (GeometryTypes::cube(dim).id() >> 1) ==
             (topologyId >> 1) &&
             (GeometryTypes::cube(dim).id() >> 1) ==
             (coerceToId >> 1)
            )>::type
        >
    {
      typedef HCube::RefinementImp<dim, CoordType> Imp;
    };
#endif

  } // namespace RefinementImp

} // namespace Dune

#endif // DUNE_GEOMETRY_REFINEMENT_HCUBE_CC
