// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_MAPPINGPROVIDER_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_MAPPINGPROVIDER_HH

#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>

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
      DUNE_EXPORT static Mapping*
      construct ( const unsigned int topologyId, const CoordVector &coords, char *mappingStorage )
      {
        static ConstructorTable< CoordVector > construct;
        return construct[ topologyId ]( coords, mappingStorage );
      }

      DUNE_EXPORT static std::size_t mappingSize ( const unsigned int topologyId )
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



    // MappingProvider
    // ---------------

    template< class ElementMapping, unsigned int codim >
    class MappingProvider;


    template< unsigned int dim, class GeometryTraits, unsigned int codim >
    class MappingProvider< HybridMapping< dim, GeometryTraits >, codim >
    {
      typedef MappingProvider< HybridMapping< dim, GeometryTraits >, codim > This;

      dune_static_assert(dim>=codim, "Codim exceeds dimension");
    public:
      static const unsigned int dimension = dim;
      static const unsigned int codimension = codim;
      static const unsigned int mydimension = dimension - codimension;

    private:
      typedef VirtualMappingFactory< mydimension, GeometryTraits > Factory;

    public:
      // Maximal amount of memory required to store a mapping
      static const unsigned int maxMappingSize = Factory::maxMappingSize;

      typedef typename Factory::Mapping Mapping;

      template< class CoordVector >
      static Mapping*
      construct ( const unsigned int topologyId, const CoordVector &coords, char *mappingStorage )
      {
        return Factory::construct( topologyId, coords, mappingStorage );
      }

      template< class CoordVector >
      static Mapping *create ( const unsigned int topologyId, const CoordVector &coords )
      {
        char *mapping = new char[ mappingSize( topologyId ) ];
        return construct( topologyId, coords, mapping );
      }

      static std::size_t mappingSize ( const unsigned int topologyId )
      {
        return Factory::mappingSize( topologyId );
      }
    };


    template< class Topology, class GeometryTraits, unsigned int codim >
    class MappingProvider< NonHybridMapping< Topology, GeometryTraits >, codim >
    {
      typedef MappingProvider< NonHybridMapping< Topology, GeometryTraits >, codim > This;

    public:
      static const unsigned int dimension = Topology::dimension;
      static const unsigned int codimension = codim;
      static const unsigned int mydimension = dimension - codimension;

      static const bool hybrid = (mydimension != dimension) && IsHybrid< Topology >::value;

    private:
      template< bool >
      struct HybridFactory
        : public VirtualMappingFactory< mydimension, GeometryTraits >
      {};

      template< bool >
      struct NonHybridFactory
        : public NonHybridMappingFactory
          < typename SubTopology< Topology, codim, 0 >::type, GeometryTraits >
      {};

      typedef typename conditional< hybrid, HybridFactory<true>, NonHybridFactory<false> >::type Factory;

    public:
      // Maximal amount of memory required to store a mapping
      static const unsigned int maxMappingSize = Factory::maxMappingSize;

      typedef typename Factory::Mapping Mapping;

      template< class CoordVector >
      static Mapping*
      construct ( const unsigned int topologyId, const CoordVector &coords, char *mappingStorage )
      {
        return Factory::construct( topologyId, coords, mappingStorage );
      }

      template< class CoordVector >
      static Mapping *create ( const unsigned int topologyId, const CoordVector &coords )
      {
        Mapping *mapping = static_cast< Mapping * >( operator new( mappingSize( topologyId ) ) );
        construct( topologyId, coords, mapping );
        return mapping;
      }

      static std::size_t mappingSize ( const unsigned int topologyId )
      {
        return Factory::mappingSize( topologyId );
      }
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_MAPPINGPROVIDER_HH
