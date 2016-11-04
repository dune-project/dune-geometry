// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_TOPOLOGYTYPES_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_TOPOLOGYTYPES_HH

#include <cassert>
#include <string>

#include <dune/common/typetraits.hh>
#include <dune/common/unused.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    enum TopologyConstruction { pyramidConstruction = 0, prismConstruction = 1 };



    // Basic Topology Types
    // --------------------

    struct Point
    {
      static const unsigned int dimension = 0;
      static const unsigned int numCorners = 1;

      static const unsigned int id = 0;

      static std :: string name ()
      {
        return "p";
      }
    };


    template< class BaseTopology >
    struct Prism
    {
      static const unsigned int dimension = BaseTopology :: dimension + 1;
      static const unsigned int numCorners = 2 * BaseTopology :: numCorners;

      static const unsigned int id = BaseTopology::id | ((unsigned int)prismConstruction << (dimension-1));

      static std :: string name ()
      {
        return BaseTopology :: name() + "l";
        // return BaseTopology :: name() + "'";
      }
    };


    template< class BaseTopology >
    struct Pyramid
    {
      static const unsigned int dimension = BaseTopology :: dimension + 1;
      static const unsigned int numCorners = BaseTopology :: numCorners + 1;

      static const unsigned int id = BaseTopology::id | ((unsigned int)pyramidConstruction << (dimension-1));

      static std :: string name ()
      {
        return BaseTopology :: name() + "o";
        // return BaseTopology :: name() + "°";
      }
    };



    template< class Topology >
    struct BaseTopology;

    template< class Base >
    struct BaseTopology< Prism< Base > >
    {
      typedef Base type;
    };

    template< class Base >
    struct BaseTopology< Pyramid< Base > >
    {
      typedef Base type;
    };



    template< class Topology >
    struct IsSimplex
    {
      static const bool value = ((Topology::id >> 1) == 0);
    };

    template< class Topology >
    struct IsCube
    {
      static const bool value = ((Topology::id | 1) == (1 << Topology::dimension) - 1);
    };

    template< class Topology >
    struct IsHybrid
    {
      static const bool value
        = !(IsSimplex< Topology >::value || IsCube< Topology >::value);
    };

    template< class Topology >
    struct IsGeneralizedPrism
    {
      static const bool value = false;
    };

    template< class BaseTopology >
    struct IsGeneralizedPrism< Prism< BaseTopology > >
    {
      static const bool value
        = (IsGeneralizedPrism< BaseTopology >::value || IsSimplex< BaseTopology >::value);
    };



    // Dynamic Topology Properties
    // ---------------------------

    /** \brief obtain the number of topologies of a given dimension
     *
     *  \note Valid topology ids are 0,...,numTopologies(dim)-1.
     *
     *  \param[in]  dim  dimension
     *
     *  \returns number of topologies for the dimension
     */
    inline static unsigned int numTopologies ( int dim ) noexcept
    {
      return (1u << dim);
    }

    /** \brief check whether a pyramid construction was used to create a given
     *         codimension
     *
     *  \param[in]  topologyId  id of the topology
     *  \param[in]  dim         dimension of the topology
     *  \param[in]  codim       codimension for which the information is desired
     *                          (defaults to 0)
     *
     *  \returns true, if a pyramid construction was used to generate the
     *           codimension the topology.
     */
    inline bool static isPyramid ( unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim < dim) );
      return (((topologyId & ~1) & (1u << (dim-codim-1))) == 0);
    }

    /** \brief check whether a prism construction was used to create a given
     *         codimension
     *
     *  \param[in]  topologyId  id of the topology
     *  \param[in]  dim         dimension of the topology
     *  \param[in]  codim       codimension for which the information is desired
     *                          (defaults to 0)
     *
     *  \returns true, if a prism construction was used to generate the
     *           codimension the topology.
     */
    inline static bool isPrism ( unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim < dim) );
      return (( (topologyId | 1) & (1u << (dim-codim-1))) != 0);
    }

    /** \brief check whether a specific topology construction was used to create a
     *         given codimension
     *
     *  \param[in]  construction  construction to check for
     *  \param[in]  topologyId    id of the topology
     *  \param[in]  dim           dimension of the topology
     *  \param[in]  codim         codimension for which the information is desired
     *                            (defaults to 0)
     *
     *  \returns true, if construction was used to generate the codimension the
     *           topology.
     */
    inline static bool isTopology ( TopologyConstruction construction, unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim <= dim) );
      return (codim >= (dim-1)) || (((topologyId >> (dim-codim-1)) & 1) == (unsigned int)construction);
    }

    /** \brief obtain the base topology of a given codimension
     *
     *  \param[in]  topologyId    id of the topology
     *  \param[in]  dim           dimension of the topology
     *  \param[in]  codim         codimension for which the information is desired
     *                            (defaults to 1)
     */
    inline static unsigned int baseTopologyId ( unsigned int topologyId, int dim, int codim = 1 ) noexcept
    {
      assert( (dim >= 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim <= dim) );
      return topologyId & ((1u << (dim-codim)) - 1);
    }



    // SimplexTopology
    // ---------------

    template< unsigned int dim >
    struct SimplexTopology
    {
      typedef Pyramid< typename SimplexTopology< dim-1 >::type > type;
    };

    template<>
    struct SimplexTopology< 0 >
    {
      typedef Point type;
    };



    // CubeTopology
    // ------------

    template< unsigned int dim >
    struct CubeTopology
    {
      typedef Prism< typename CubeTopology< dim-1 >::type > type;
    };

    template<>
    struct CubeTopology< 0 >
    {
      typedef Point type;
    };



    // PyramidTopology
    // --------------

    template< unsigned int dim >
    struct PyramidTopology
    {
      typedef Pyramid< typename CubeTopology< dim-1 >::type > type;
    };



    // PrismTopology
    // --------------

    template< unsigned int dim >
    struct PrismTopology
    {
      typedef Prism< typename SimplexTopology< dim-1 >::type > type;
    };




    // IfTopology
    // ----------

    template< template< class > class Operation, int dim, class Topology = Point >
    class IfTopology
    {
      typedef IfTopology< Operation, dim-1, Prism< Topology > > IfPrism;
      typedef IfTopology< Operation, dim-1, Pyramid< Topology > > IfPyramid;

    public:
      template< class... Args >
      static void apply ( unsigned int topologyId, Args &&... args )
      {
        if( topologyId & 1 )
          IfPrism::apply( topologyId >> 1, std::forward< Args >( args )... );
        else
          IfPyramid::apply( topologyId >> 1, std::forward< Args >( args )... );
      }
    };

    template< template< class > class Operation, class Topology >
    class IfTopology< Operation, 0, Topology >
    {
    public:
      template< class... Args >
      static void apply ( unsigned int topologyId, Args &&... args )
      {
        DUNE_UNUSED_PARAMETER( topologyId );
        Operation< Topology >::apply( std::forward< Args >( args )... );
      }
    };

  }

}

#endif // DUNE_GEOMETRY_GENERICGEOMETRY_TOPOLOGYTYPES_HH
