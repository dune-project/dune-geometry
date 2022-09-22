// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DEPRECATED_TOPOLOGY_HH
#define DUNE_DEPRECATED_TOPOLOGY_HH

  namespace Impl
  {

    // Basic Topology Types
    // --------------------

    // PointDeprecationHelper can be used to prevent a deprecation warning for Point
    struct PointDeprecationHelper
    {
      static const unsigned int dimension = 0;
      static const unsigned int numCorners = 1;

      static const unsigned int id = 0;

      static std::string name () { return "p"; }
    };

    using Point [[deprecated("Use GeometryTypes::vertex instead.")]] = PointDeprecationHelper;


    template< class BaseTopology >
    struct [[deprecated("Use GeometryTypes::prismaticExtension(GeometryType gt) instead.")]] Prism
    {
      static const unsigned int dimension = BaseTopology::dimension + 1;
      static const unsigned int numCorners = 2 * BaseTopology::numCorners;

      static const unsigned int id = BaseTopology::id | ((unsigned int)prismConstruction << (dimension-1));

      static std::string name () { return BaseTopology::name() + "l"; }
    };


    template< class BaseTopology >
    struct [[deprecated("Use GeometryTypes::conicalExtension(GeometryType gt) instead.")]] Pyramid
    {
      static const unsigned int dimension = BaseTopology::dimension + 1;
      static const unsigned int numCorners = BaseTopology::numCorners + 1;

      static const unsigned int id = BaseTopology::id | ((unsigned int)pyramidConstruction << (dimension-1));

      static std::string name () { return BaseTopology::name() + "o"; }
    };



    // Properties of Topologies
    // ------------------------

    template< class Topology >
    struct [[deprecated("Use GeometryType::isSimplex() instead.")]] IsSimplex
      : public std::integral_constant< bool, (Topology::id >> 1) == 0 >
    {};

    template< class Topology >
    struct [[deprecated("Use GeometryType::isCube() instead.")]] IsCube
      : public std::integral_constant< bool,  (Topology::id | 1) == (1 << Topology::dimension) - 1 >
    {};

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
    [[deprecated("Use GeometryType::isPrismatic() or GeometryType::isConical() instead.")]]
    inline static bool isTopology ( TopologyConstruction construction, unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim <= dim) );
      return (codim >= (dim-1)) || (((topologyId >> (dim-codim-1)) & 1) == (unsigned int)construction);
    }


    // SimplexTopology
    // ---------------

    template< unsigned int dim >
    struct [[deprecated("Use GeometryTypes::simplex(dim) instead.")]] SimplexTopology
    {
      typedef Pyramid< typename SimplexTopology< dim-1 >::type > type;
    };

    template<>
    struct [[deprecated("Use GeometryTypes::simplex(dim) instead.")]] SimplexTopology< 0 >
    {
      typedef Point type;
    };



    // CubeTopology
    // ------------

    template< unsigned int dim >
    struct [[deprecated("Use GeometryTypes::cube(dim) instead.")]] CubeTopology
    {
      typedef Prism< typename CubeTopology< dim-1 >::type > type;
    };

    template<>
    struct [[deprecated("Use GeometryTypes::simplex(dim) instead.")]] CubeTopology< 0 >
    {
      typedef Point type;
    };



    // PyramidTopology
    // ---------------

    template< unsigned int dim >
    struct [[deprecated]] PyramidTopology
    {
      typedef Pyramid< typename CubeTopology< dim-1 >::type > type;
    };



    // PrismTopology
    // -------------

    template< unsigned int dim >
    struct [[deprecated]] PrismTopology
    {
      typedef Prism< typename SimplexTopology< dim-1 >::type > type;
    };




    // IfTopology
    // ----------

    template< template< class > class Operation, int dim, class Topology = PointDeprecationHelper >
    struct [[deprecated("Use IfGeometryType instead.")]] IfTopology
    {
      template< class... Args >
      static auto apply ( unsigned int topologyId, Args &&... args )
      {
        if( topologyId & 1 )
          return IfTopology< Operation, dim-1, Prism< Topology > >::apply( topologyId >> 1, std::forward< Args >( args )... );
        else
          return IfTopology< Operation, dim-1, Pyramid< Topology > >::apply( topologyId >> 1, std::forward< Args >( args )... );
      }
    };

    template< template< class > class Operation, class Topology >
    struct [[deprecated("Use IfGeometryType instead.")]] IfTopology< Operation, 0, Topology >
    {
      template< class... Args >
      static auto apply ([[maybe_unused]] unsigned int topologyId, Args &&... args)
      {
        return Operation< Topology >::apply( std::forward< Args >( args )... );
      }
    };

  } // namespace Impl
#endif
