// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_REFINEMENT_HCUBETRIANGULATION_CC
#define DUNE_GEOMETRY_REFINEMENT_HCUBETRIANGULATION_CC

/*!
 * \file
 * \brief This file contains the \ref Refinement implementation for
 *        triangulating hypercubes (quadrilateral -> triangle,
 *        hexahedron -> tetrahedron)
 *
 * See \ref HCubeTriangulation.
 */

/*!
 * \defgroup HCubeTriangulation Refinement implementation for triangulating hypercubes
 * \ingroup Refinement
 *
 * Most stuff here is explained in \ref SimplexRefinement.
 *
 * We simply triangulate the hypercube into Kuhn simplices, and than
 * use the \ref SimplexRefinement to do the refinement.
 *
 * We explicitly use some of the utilities from the \ref SimplexRefinement.
 */

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include "base.cc"
#include "simplex.cc"

namespace Dune
{
  namespace RefinementImp
  {
    /*!
     * \brief This namespace contains the \ref Refinement
     *        implementation for triangulating hypercubes
     *        (GeometryType::cube -> GeometryType::simplex)
     *
     * See \ref HCubeTriangulation.
     */
    namespace HCubeTriangulation {

      // ////////////
      //
      //  Utilities
      //

      using Simplex::getPermutation;
      using Simplex::referenceToKuhn;

      // ////////////////////////////////////
      //
      //  Refine a hypercube with simplices
      //

      // forward declaration of the iterator base
      template<int dimension, class CoordType, int codimension>
      class RefinementIteratorSpecial;

      template<int dimension_, class CoordType>
      class RefinementImp
      {
      public:
        constexpr static int dimension = dimension_;

        typedef CoordType ctype;

        template<int codimension>
        struct Codim;
        typedef typename Codim<dimension>::SubEntityIterator VertexIterator;
        typedef FieldVector<CoordType, dimension> CoordVector;
        typedef typename Codim<0>::SubEntityIterator ElementIterator;
        typedef FieldVector<int, dimension+1> IndexVector;

        static int nVertices(int nIntervals);
        static VertexIterator vBegin(int nIntervals);
        static VertexIterator vEnd(int nIntervals);

        static int nElements(int nIntervals);
        static ElementIterator eBegin(int nIntervals);
        static ElementIterator eEnd(int nIntervals);
      private:
        friend class RefinementIteratorSpecial<dimension, CoordType, 0>;
        friend class RefinementIteratorSpecial<dimension, CoordType, dimension>;

        typedef Simplex::RefinementImp<dimension, CoordType> BackendRefinement;
      };

      template<int dimension, class CoordType>
      template<int codimension>
      struct RefinementImp<dimension, CoordType>::Codim
      {
        class SubEntityIterator;
        typedef Dune::MultiLinearGeometry<CoordType,dimension-codimension,dimension> Geometry;
      };

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nVertices(int nIntervals)
      {
        return BackendRefinement::nVertices(nIntervals) * factorial(int(dimension));
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vBegin(int nIntervals)
      {
        return VertexIterator(nIntervals);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vEnd(int nIntervals)
      {
        return VertexIterator(nIntervals, true);
      }

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nElements(int nIntervals)
      {
        return BackendRefinement::nElements(nIntervals) * factorial(int(dimension));
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eBegin(int nIntervals)
      {
        return ElementIterator(nIntervals);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eEnd(int nIntervals)
      {
        return ElementIterator(nIntervals, true);
      }

      // //////////////
      //
      // The iterator
      //

      // vertices
      template<int dimension, class CoordType>
      class RefinementIteratorSpecial<dimension, CoordType, dimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::CoordVector CoordVector;
        typedef typename Refinement::template Codim<dimension>::Geometry Geometry;

        RefinementIteratorSpecial(int nIntervals, bool end = false);

        void increment();

        CoordVector coords() const;

        Geometry geometry() const;

        int index() const;
      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<dimension>::SubEntityIterator BackendIterator;
        constexpr static int nKuhnSimplices = factorial(int(dimension));

        int nIntervals_;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      RefinementIteratorSpecial(int nIntervals, bool end)
        : nIntervals_(nIntervals), kuhnIndex(0),
          backend(BackendRefinement::vBegin(nIntervals_)),
          backendEnd(BackendRefinement::vEnd(nIntervals_))
      {
        if (end)
          kuhnIndex = nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      increment()
      {
        ++backend;
        if (backend == backendEnd)
        {
          backend = BackendRefinement::vBegin(nIntervals_);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      coords() const
      {
        return referenceToKuhn(backend.coords(), getPermutation<dimension>(kuhnIndex));
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::Geometry
      RefinementIteratorSpecial<dimension, CoordType, dimension>::geometry () const
      {
        std::vector<CoordVector> corners(1);
        corners[0] = referenceToKuhn(backend.coords(), getPermutation<dimension>(kuhnIndex));
        return Geometry(GeometryTypes::vertex, corners);
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nVertices(nIntervals_) + backend.index();
      }

      // elements
      template<int dimension, class CoordType>
      class RefinementIteratorSpecial<dimension, CoordType, 0>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::IndexVector IndexVector;
        typedef typename Refinement::CoordVector CoordVector;
        typedef typename Refinement::template Codim<0>::Geometry Geometry;

        RefinementIteratorSpecial(int nIntervals_, bool end = false);
        RefinementIteratorSpecial(const RefinementIteratorSpecial<dimension, CoordType, 0> &other);

        void increment();

        IndexVector vertexIndices() const;
        int index() const;
        CoordVector coords() const;

        Geometry geometry() const;

      private:
        CoordVector global(const CoordVector &local) const;

      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<0>::SubEntityIterator BackendIterator;
        constexpr static int nKuhnSimplices = factorial(dimension);

        int nIntervals_;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      RefinementIteratorSpecial(int nIntervals, bool end)
        : nIntervals_(nIntervals), kuhnIndex(0),
          backend(BackendRefinement::eBegin(nIntervals_)),
          backendEnd(BackendRefinement::eEnd(nIntervals_))
      {
        if (end)
          kuhnIndex = nKuhnSimplices;
      }
      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      RefinementIteratorSpecial(const RefinementIteratorSpecial<dimension, CoordType, 0> &other)
        : nIntervals_(other.nIntervals_), kuhnIndex(other.kuhnIndex),
          backend(other.backend),
          backendEnd(other.backendEnd)
      {}

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      increment()
      {
        ++backend;
        if (backend == backendEnd)
        {
          backend = BackendRefinement::eBegin(nIntervals_);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::IndexVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      vertexIndices() const
      {
        IndexVector indices = backend.vertexIndices();

        int base = kuhnIndex * BackendRefinement::nVertices(nIntervals_);
        indices += base;

        return indices;
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nElements(nIntervals_) + backend.index();
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      coords() const
      {
        return global(backend.coords());
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::Geometry
      RefinementIteratorSpecial<dimension, CoordType, 0>::geometry () const
      {
        const typename BackendIterator::Geometry &bgeo =
          backend.geometry();
        std::vector<CoordVector> corners(dimension+1);
        for(int i = 0; i <= dimension; ++i)
          corners[i] = global(bgeo.corner(i));

        return Geometry(bgeo.type(), corners);
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      global(const CoordVector &local) const
      {
        return referenceToKuhn(local, getPermutation<dimension>(kuhnIndex));
      }

      // common
      template<int dimension, class CoordType>
      template<int codimension>
      class RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator
        : public ForwardIteratorFacade<typename RefinementImp<dimension, CoordType>::template Codim<codimension>::SubEntityIterator, int>,
          public RefinementIteratorSpecial<dimension, CoordType, codimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef SubEntityIterator This;

        SubEntityIterator(int nIntervals, bool end = false);

        bool equals(const This &other) const;
      protected:
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::kuhnIndex;
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::backend;
      };

#ifndef DOXYGEN
      template<int dimension, class CoordType>
      template<int codimension>
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      SubEntityIterator(int nIntervals, bool end)
        : RefinementIteratorSpecial<dimension, CoordType, codimension>(nIntervals, end)
      {}

      template<int dimension, class CoordType>
      template<int codimension>
      bool
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      equals(const This &other) const
      { return kuhnIndex == other.kuhnIndex && backend == other.backend; }

#endif // DOXYGEN

    } // namespace HCubeTriangulation
  } // namespace RefinementImp

  namespace RefinementImp
  {
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
             (GeometryTypes::simplex(dim).id() >> 1) ==
             (coerceToId >> 1)
            )>::type
        >
    {
      typedef HCubeTriangulation::RefinementImp<dim, CoordType> Imp;
    };
#endif

  } // namespace RefinementImp
} // namespace Dune

#endif // DUNE_GEOMETRY_REFINEMENT_HCUBETRIANGULATION_CC
