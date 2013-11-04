// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_CONVERSION_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_CONVERSION_HH

#include <dune/common/static_assert.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // DuneGeometryType
    // ----------------

    /** \class DuneGeometryType
     *  \ingroup GenericGeometry
     *  \brief statically convert a generic topology type into a GeometryType
     *
     *  \tparam  Topology  topology type to be converted
     *  \tparam  linetype  basic geometry type to assign to a line
     *                     (either simplex or cube)
     */
    template< class Topology, GeometryType::BasicType linetype >
    class DuneGeometryType;

    template< GeometryType::BasicType linetype >
    class DuneGeometryType< Point, linetype >
    {
      dune_static_assert( (linetype == GeometryType::simplex)
                          || (linetype == GeometryType::cube),
                          "Parameter linetype may only be a simplex or a cube." );

    public:
      static const unsigned int dimension = 0;
      static const GeometryType::BasicType basicType = linetype;
    };

    template< class BaseTopology, GeometryType::BasicType linetype >
    class DuneGeometryType< Prism< BaseTopology >, linetype >
    {
      typedef DuneGeometryType< BaseTopology, linetype > DuneBaseGeometryType;

      dune_static_assert( (linetype == GeometryType::simplex)
                          || (linetype == GeometryType::cube),
                          "Parameter linetype may only be a simplex or a cube." );

      dune_static_assert( (DuneBaseGeometryType::basicType == GeometryType::simplex)
                          || (DuneBaseGeometryType::basicType == GeometryType::cube),
                          "Only prisms over simplices or cubes can be converted." );

    public:
      static const unsigned int dimension = DuneBaseGeometryType::dimension + 1;
      static const GeometryType::BasicType basicType
        = ((dimension == 1)
           ? linetype
           : ((dimension == 2) || (DuneBaseGeometryType::basicType == GeometryType::cube))
           ? GeometryType::cube
           : GeometryType::prism);
    };

    template< class BaseTopology, GeometryType::BasicType linetype >
    class DuneGeometryType< Pyramid< BaseTopology >, linetype >
    {
      typedef DuneGeometryType< BaseTopology, linetype > DuneBaseGeometryType;

      dune_static_assert( (linetype == GeometryType::simplex)
                          || (linetype == GeometryType::cube),
                          "Parameter linetype may only be a simplex or a cube." );

      dune_static_assert( (DuneBaseGeometryType::basicType == GeometryType::simplex)
                          || (DuneBaseGeometryType::basicType == GeometryType::cube),
                          "Only pyramids over simplices or cubes can be converted." );

    public:
      static const unsigned int dimension = DuneBaseGeometryType::dimension + 1;
      static const GeometryType::BasicType basicType
        = ((dimension == 1)
           ? linetype
           : ((dimension == 2) || (DuneBaseGeometryType::basicType == GeometryType::simplex))
           ? GeometryType::simplex
           : GeometryType::pyramid);
    };



    // DuneGeometryTypeProvider
    // ------------------------

    /** \class DuneGeometryTypeProvider
     *  \brief dynamically convert a generic topology type into a GeometryType
     *
     *  \tparam  dim       dimension of the topologies to be converted
     *  \tparam  linetype  basic geometry type to assign to a line
     *                     (either simplex or cube)
     *
     *  \note After 3D not all geometries are simplices, pyramids, prisms or
     *        cubes so that no meaningful GeometryType is available; therefore
     *        none is returned.
     */
    template< unsigned int dim, GeometryType::BasicType linetype >
    struct DuneGeometryTypeProvider
    {
      /** \brief dimension of the topologies to be converted */
      static const unsigned int dimension = dim;

      /** \brief number of possible topologies */
      static const unsigned int numTopologies = (1 << dimension);

    private:
      GeometryType types_[ (dimension>=1) ? numTopologies / 2 : numTopologies ];

      DUNE_EXPORT static const DuneGeometryTypeProvider &instance ()
      {
        static DuneGeometryTypeProvider inst;
        return inst;
      }

    public:
      /** \brief obtain a Geometry type from a topology id
       *
       *  \param[in]  topologyId  id of the topology to be converted
       *
       *  \returns GeometryType associated with the given topology id
       */
      static const GeometryType &type ( unsigned int topologyId )
      {
        assert( topologyId < numTopologies );
        return instance().types_[ topologyId / 2 ];
      }
    };


    // Convert
    // -------

    template< GeometryType :: BasicType type, unsigned int dim >
    struct Convert;

    template< unsigned int dim >
    struct Convert< GeometryType :: simplex, dim >
    {
      typedef Pyramid
      < typename Convert< GeometryType :: simplex, dim-1 > :: type >
      type;
    };

    template<>
    struct Convert< GeometryType :: simplex, 0 >
    {
      typedef Point type;
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: cube, dim >
    {
      typedef Prism< typename Convert< GeometryType :: cube, dim-1 > :: type >
      type;
    };

    template<>
    struct Convert< GeometryType :: cube, 0 >
    {
      typedef Point type;
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: prism, dim >
    {
      typedef Prism
      < typename Convert< GeometryType :: simplex, dim-1 > :: type >
      type;
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: pyramid, dim >
    {
      typedef Pyramid
      < typename Convert< GeometryType :: cube, dim-1 > :: type >
      type;
    };

  }

}

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_CONVERSION_HH
