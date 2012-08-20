// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_MAPPINGPROVIDER_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_MAPPINGPROVIDER_HH

#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/maximum.hh>
#include <dune/geometry/genericgeometry/cachedmapping.hh>
#include <dune/geometry/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // NonHybridMappingFactory
    // -----------------------

    template< class Topology, class GeometryTraits >
    class NonHybridMappingFactory
    {
      typedef NonHybridMappingFactory< Topology, GeometryTraits > This;

    public:
      typedef NonHybridMapping< Topology, GeometryTraits > Mapping;

      static const unsigned int maxMappingSize = sizeof( Mapping );

      template< class CoordVector >
      static Mapping*
      construct ( const unsigned int topologyId, const CoordVector &coords, char *mappingStorage )
      {
        assert( (topologyId >> 1) == (Topology::id >> 1) );
        return new( mappingStorage ) Mapping( coords );
      }

      static std::size_t mappingSize ( const unsigned int topologyId )
      {
        return sizeof( Mapping );
      }
    };



    // VirtualMappingFactory
    // ---------------------

    template< unsigned int dim, class GeometryTraits >
    class VirtualMappingFactory
    {
      typedef VirtualMappingFactory< dim, GeometryTraits > This;

      static const unsigned int numTopologies = (1 << dim);

      template< int topologyId >
      struct MappingSize
      {
        typedef typename GenericGeometry::Topology< (unsigned int) topologyId, dim >::type Topology;
        static const int v = sizeof( VirtualMapping< Topology, GeometryTraits > );

        static void apply ( std::size_t (&mappingSize)[ numTopologies ] )
        {
          mappingSize[ topologyId ] = v;
        }
      };

      template< class CoordVector >
      class ConstructorTable;

      struct MappingSizeCache;

    public:
      typedef HybridMapping< dim, GeometryTraits > Mapping;

      static const unsigned int maxMappingSize = Maximum< MappingSize, 0, ((numTopologies > 0) ? (numTopologies-1) : 0) >::v;

      template< class CoordVector >
      static Mapping*
      construct ( const unsigned int topologyId, const CoordVector &coords, char *mappingStorage )
      {
        static ConstructorTable< CoordVector > construct;
        return construct[ topologyId ]( coords, mappingStorage );
      }

      static std::size_t mappingSize ( const unsigned int topologyId )
      {
        static MappingSizeCache mappingSize;
        return mappingSize[ topologyId ];
      }
    };


    // VirtualMappingFactory::ConstructorTable
    // ---------------------------------------

    template< unsigned int dim, class GeometryTraits >
    template< class CoordVector >
    class VirtualMappingFactory< dim, GeometryTraits >::ConstructorTable
    {
      typedef Mapping* (*Construct)( const CoordVector &coords, char *mappingStorage );

      template< int i >
      struct Builder;

    public:
      ConstructorTable ()
      {
        ForLoop< Builder, 0, numTopologies-1 >::apply( construct_ );
      }

      Construct operator[] ( const unsigned int topologyId )
      {
        assert( topologyId < numTopologies );
        return construct_[ topologyId ];
      }

    private:
      template< class Topology >
      static Mapping*
      construct ( const CoordVector &coords, char *mappingStorage )
      {
        typedef VirtualMapping< Topology, GeometryTraits > VMapping;
        return new( mappingStorage ) VMapping( coords );
      }

      Construct construct_[ numTopologies ];
    };


    // VirtualMappingFactory::ConstructorTable::Builder
    // ------------------------------------------------

    template< unsigned int dim, class GeometryTraits >
    template< class CoordVector >
    template< int topologyId >
    struct VirtualMappingFactory< dim, GeometryTraits >::ConstructorTable< CoordVector >::Builder
    {
      static void apply ( Construct (&construct)[ numTopologies ] )
      {
        typedef typename GenericGeometry::Topology< (unsigned int) topologyId, dim >::type Topology;
        construct[ topologyId ] = ConstructorTable< CoordVector >::template construct< Topology >;
      }
    };



    // VirtualMappingFactory::MappingSizeCache
    // ---------------------------------------

    template< unsigned int dim, class GeometryTraits >
    struct VirtualMappingFactory< dim, GeometryTraits >::MappingSizeCache
    {
      MappingSizeCache ()
      {
        ForLoop< MappingSize, 0, numTopologies-1 >::apply( size_ );
      }

      std::size_t operator[] ( const unsigned int topologyId )
      {
        assert( topologyId < numTopologies );
        return size_[ topologyId ];
      }

    private:
      std::size_t size_[ numTopologies ];
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_MAPPINGPROVIDER_HH
