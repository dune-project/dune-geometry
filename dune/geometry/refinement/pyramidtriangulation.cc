// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_REFINEMENT_PYRAMIDTRIANGULATION_CC
#define DUNE_GEOMETRY_REFINEMENT_PYRAMIDTRIANGULATION_CC

#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

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
     * implementation for triangulating pyramids
     * (GeometryType::pyramid -> GeometryType::simplex)
     */
    namespace PyramidTriangulation
    {
      // ////////////
      //
      //  Utilities
      //

      using Simplex::getPermutation;
      using Simplex::referenceToKuhn;

      // ////////////////////////////////////
      //
      //  Refine a pyramid with simplices
      //

      // forward declaration of the iterator base
      template<int dimension, class CoordType, int codimension>
      class RefinementIteratorSpecial;

      /*
       * The permutations 0 and 1 of the Kuhn-decomposition of a cube into simplices form a pyramid.
       * The resulting pyramid is not oriented the same as the reference pyramid and so the Kuhn-coordinates
       * have to be transformed using the method below.
       */
      template<int dimension, class CoordType> FieldVector<CoordType, dimension>
      transformCoordinate(/*! Point to transform */ FieldVector<CoordType, dimension> point)
      {
        FieldVector<CoordType, dimension> transform;
        transform[0]=1-point[0];
        transform[1]=1-point[1];
        transform[2]=point[2];
        return transform;
      }

      /*!
       * \brief Implementation of the refinement of a pyramid into simplices.
       *
       * Note that the virtual vertices of two intersecting simplices might have copies, i.e.
       * by running over all vertices using the VertexIterator you might run over some twice.
       */
      template<int dimension_, class CoordType>
      class RefinementImp
      {
      public:
        enum {dimension = dimension_};

        typedef CoordType ctype;

        template<int codimension>
        struct Codim;
        typedef typename Codim<dimension>::SubEntityIterator VertexIterator;
        typedef FieldVector<CoordType, dimension> CoordVector;
        typedef typename Codim<0>::SubEntityIterator ElementIterator;
        typedef FieldVector<int, dimension+1> IndexVector;

        static int nVertices(int nPyramids);
        static VertexIterator vBegin(int nPyramids);
        static VertexIterator vEnd(int nPyramids);

        static int nElements(int nPyramids);
        static ElementIterator eBegin(int nPyramids);
        static ElementIterator eEnd(int nPyramids);

      private:
        friend class RefinementIteratorSpecial<dimension, CoordType, 0>;
        friend class RefinementIteratorSpecial<dimension, CoordType, dimension>;

        typedef Simplex::RefinementImp<dimension, CoordType> BackendRefinement;

        enum { nKuhnSimplices = 2 };
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
      nVertices(int nPyramids)
      {
        return BackendRefinement::nVertices(nPyramids) * nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vBegin(int nPyramids)
      {
        return VertexIterator(nPyramids);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vEnd(int nPyramids)
      {
        return VertexIterator(nPyramids, true);
      }

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nElements(int nPyramids)
      {
        return BackendRefinement::nElements(nPyramids) * nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eBegin(int nPyramids)
      {
        return ElementIterator(nPyramids);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eEnd(int nPyramids)
      {
        return ElementIterator(nPyramids, true);
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

        RefinementIteratorSpecial(int nPyramids, bool end = false);

        void increment();

        CoordVector coords() const;

        Geometry geometry() const;

        int index() const;
      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<dimension>::SubEntityIterator BackendIterator;
        enum { nKuhnSimplices = 2 };

        int nPyramids_;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      RefinementIteratorSpecial(int nPyramids, bool end)
        : nPyramids_(nPyramids), kuhnIndex(0),
          backend(BackendRefinement::vBegin(nPyramids_)),
          backendEnd(BackendRefinement::vEnd(nPyramids_))
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
        if(backend == backendEnd)
        {
          backend = BackendRefinement::vBegin(nPyramids_);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      coords() const
      {
        return transformCoordinate(referenceToKuhn(backend.coords(),
                                                   getPermutation<dimension>(kuhnIndex)));
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::Geometry
      RefinementIteratorSpecial<dimension, CoordType, dimension>::geometry () const
      {
        std::vector<CoordVector> corners(1);
        corners[0] = referenceToKuhn(backend.coords(), getPermutation<dimension>(kuhnIndex));
        return Geometry(GeometryType(0), corners);
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nVertices(nPyramids_) + backend.index();
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

        RefinementIteratorSpecial(int nPyramids, bool end = false);

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
        enum { nKuhnSimplices = 2};

        int nPyramids_;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      RefinementIteratorSpecial(int nPyramids, bool end)
        : nPyramids_(nPyramids), kuhnIndex(0),
          backend(BackendRefinement::eBegin(nPyramids_)),
          backendEnd(BackendRefinement::eEnd(nPyramids_))
      {
        if (end)
          kuhnIndex = nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      increment()
      {
        ++backend;
        if (backend == backendEnd)
        {
          backend = BackendRefinement::eBegin(nPyramids_);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::IndexVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      vertexIndices() const
      {
        IndexVector indices = backend.vertexIndices();

        int base = kuhnIndex * BackendRefinement::nVertices(nPyramids_);
        indices += base;

        return indices;
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nElements(nPyramids_) + backend.index();
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
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      geometry() const
      {
        const typename BackendIterator::Geometry &
        bgeo = backend.geometry();
        std::vector<CoordVector> corners(dimension+1);
        for(int i = 0; i <= dimension; ++i)
          corners[i] = global(bgeo.corner(i));

        return Geometry(bgeo.type(), corners);
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::
      CoordVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      global(const CoordVector &local) const
      {
        return transformCoordinate(referenceToKuhn(local,
                                                   getPermutation<dimension>(kuhnIndex)));
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

        SubEntityIterator(int nPyramids, bool end = false);

        bool equals(const This &other) const;
      protected:
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::kuhnIndex;
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::backend;
      };

#ifndef DOXYGEN
      template<int dimension, class CoordType>
      template<int codimension>
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      SubEntityIterator(int nPyramids, bool end)
        : RefinementIteratorSpecial<dimension, CoordType, codimension>(nPyramids, end)
      {}

      template<int dimension, class CoordType>
      template<int codimension>
      bool
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      equals(const This &other) const
      {
        return kuhnIndex == other.kuhnIndex && backend == other.backend;
      }
#endif

    } // namespace PyramidTriangulation
  } // namespace RefinementImp

  namespace RefinementImp
  {
    // ///////////////////////
    //
    // The refinement traits
    //
#ifndef DOXYGEN
    template<unsigned topologyId, class CoordType, unsigned coerceToId>
    struct Traits<
        topologyId, CoordType, coerceToId, 3,
        typename std::enable_if<
            (Impl::PyramidTopology<3>::type::id >> 1) ==
            (topologyId >> 1) &&
            (Impl::SimplexTopology<3>::type::id >> 1) ==
            (coerceToId >> 1)
            >::type>
    {
      typedef PyramidTriangulation::RefinementImp<3, CoordType> Imp;
    };
#endif

  } // namespace RefinementImp
} // namespace Dune

#endif // DUNE_GEOMETRY_REFINEMENT_PYRAMIDTRIANGULATION_CC
