// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTS_HH

#include <dune/common/forloop.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/subtopologies.hh>
#include <dune/geometry/genericgeometry/referencedomain.hh>
#include <dune/geometry/genericgeometry/hybridmapping.hh>
#include <dune/geometry/genericgeometry/mappingprovider.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElement;

  template< class ctype, int dim >
  class ReferenceElementContainer;



  // ReferenceElement for ctype = void
  // ---------------------------------

  template< int dim >
  class ReferenceElement< void, dim >
  {
    typedef ReferenceElement< void, dim > This;

    friend class ReferenceElementContainer< void, dim >;

    struct SubEntityInfo;

    // make copy constructor private
    ReferenceElement ( const This & );

  protected:
    ReferenceElement () {}
    ~ReferenceElement () {}

  public:
    /** \brief number of subentities of codimension c
     *
     *  \param[in]  c  codimension whose size is desired
     */
    int size ( int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ].size();
    }

    /** \brief number of subentities of codimension cc of subentity (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. This method returns the number of subentities
     *  of codimension cc of the current reference element, that are also
     *  a subentity of E.
     *
     *  \param[in]  i   number of subentity E (0 <= i < size( c ))
     *  \param[in]  c   codimension of subentity E
     *  \param[in]  cc  codimension whose size is desired (c <= cc <= dim)
     */
    int size ( int i, int c, int cc ) const
    {
      assert( (i >= 0) && (i < size( c )) );
      return info_[ c ][ i ].size( cc );
    }

    /** \brief obtain number of ii-th subentity with codim cc of (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. And denote by S the ii-th subentity of codimension
     *  (cc-c) of E. Then, S is a also a subentity of codimension c of the current
     *  reference element. This method returns the number of S with respect
     *  to the current reference element.
     *
     *  \param[in]  i   number of subentity E (0 <= i < size( c ))
     *  \param[in]  c   codimension of subentity E
     *  \param[in]  ii  number of subentity S (with respect to E)
     *  \param[in]  cc  codimension of subentity S (c <= cc <= dim)
     */
    int subEntity ( int i, int c, int ii, int cc ) const
    {
      assert( (i >= 0) && (i < size( c )) );
      return info_[ c ][ i ].number( ii, cc );
    }

    /** \brief obtain the type of subentity (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. This method returns the GeometryType of E.
     *
     *  \param[in]  i      number of subentity E (0 <= i < size( c ))
     *  \param[in]  c      codimension of subentity E
     */
    const GeometryType &type ( int i, int c ) const
    {
      assert( (i >= 0) && (i < size( c )) );
      return info_[ c ][ i ].type();
    }

    /** \brief obtain the type of this reference element */
    const GeometryType &type () const { return type( 0, 0 ); }

    /** \brief initialize the reference element
     *
     *  \param[in]  topologyId  topology id for the desired reference element
     */
    void initializeTopology ( unsigned int topologyId )
    {
      assert( topologyId < GenericGeometry::numTopologies( dim ) );

      // set up subentities
      for( int codim = 0; codim <= dim; ++codim )
      {
        const unsigned int size = GenericGeometry::size( topologyId, dim, codim );
        info_[ codim ].resize( size );
        for( unsigned int i = 0; i < size; ++i )
          info_[ codim ][ i ].initialize( topologyId, codim, i );
      }
    }

    /** \brief initialize the reference element
     *
     *  \tparam  Topology  topology of the desired reference element
     *
     *  \note The dimension of the topology must match dim.
     */
    template< class Topology >
    void initializeTopology ()
    {
      dune_static_assert( (Topology::dimension == dim),
                          "Cannot initialize reference element for different dimension." );
      initializeTopology( Topology::id );
    }

  private:
    std::vector< SubEntityInfo > info_[ dim+1 ];
  };



  // ReferenceElement::SubEntityInfo
  // -------------------------------

  /** \brief topological information about the subentities of a reference element */
  template< int dim >
  struct ReferenceElement< void, dim >::SubEntityInfo
  {
    SubEntityInfo ()
      : numbering_( nullptr )
    {
      std::fill( offset_, offset_ + (dim+2), 0 );
    }

    SubEntityInfo ( const SubEntityInfo &other )
      : type_( other.type_ )
    {
      std::copy( other.offset_, other.offset_ + (dim+2), offset_ );
      numbering_ = allocate();
      std::copy( other.numbering_, other.numbering_ + capacity(), numbering_ );
    }

    ~SubEntityInfo () { deallocate( numbering_ ); }

    const SubEntityInfo &operator= ( const SubEntityInfo &other )
    {
      type_ = other.type_;
      std::copy( other.offset_, other.offset_ + (dim+2), offset_ );

      deallocate( numbering_ );
      numbering_ = allocate();
      std::copy( other.numbering_, other.numbering_ + capacity(), numbering_ );

      return *this;
    }

    int size ( int cc ) const
    {
      assert( (cc >= codim()) && (cc <= dim) );
      return (offset_[ cc+1 ] - offset_[ cc ]);
    }

    int number ( int ii, int cc ) const
    {
      assert( (ii >= 0) && (ii < size( cc )) );
      return numbering_[ offset_[ cc ] + ii ];
    }

    const GeometryType &type () const { return type_; }

    void initialize ( unsigned int topologyId, int codim, unsigned int i )
    {
      const unsigned int subId = GenericGeometry::subTopologyId( topologyId, dim, codim, i );
      type_ = GeometryType( subId, dim-codim );

      // compute offsets
      for( int cc = 0; cc <= codim; ++cc )
        offset_[ cc ] = 0;
      for( int cc = codim; cc <= dim; ++cc )
        offset_[ cc+1 ] = offset_[ cc ] + GenericGeometry::size( subId, dim-codim, cc-codim );

      // compute subnumbering
      deallocate( numbering_ );
      numbering_ = allocate();
      for( int cc = codim; cc <= dim; ++cc )
        GenericGeometry::subTopologyNumbering( topologyId, dim, codim, i, cc-codim, numbering_+offset_[ cc ], numbering_+offset_[ cc+1 ] );
    }

  protected:
    int codim () const { return dim - type().dim(); }

    unsigned int *allocate () { return (capacity() != 0 ? new unsigned int[ capacity() ] : nullptr); }
    void deallocate ( unsigned int *ptr ) { delete[] ptr; }
    unsigned int capacity () const { return offset_[ dim+1 ]; }

  private:
    unsigned int *numbering_;
    unsigned int offset_[ dim+2 ];
    GeometryType type_;
  };



  // ReferenceElement
  // ----------------

  /** \class ReferenceElement
   *  \ingroup GeometryReferenceElements
   *  \brief This class provides access to geometric and topological
   *  properties of a reference element. This includes its type,
   *  the number of subentities, the volume, and a method for checking
   *  if a point is inside.
   *  The embedding of each subentity into the reference element is also
   *  provided.
   *
   *  A singleton of this class for a given geometry type can be accessed
   *  through the ReferenceElements class.

   *  \tparam ctype  field type for coordinates
   *  \tparam dim    dimension of the reference element
   *
   */
  template< class ctype, int dim >
  class ReferenceElement
    : public ReferenceElement< void, dim >
  {
    typedef ReferenceElement< ctype, dim > This;
    typedef ReferenceElement< void, dim > Base;

    friend class ReferenceElementContainer< ctype, dim >;

    // make copy constructor private
    ReferenceElement ( const This & );

    ReferenceElement () {}

    ~ReferenceElement ()
    {
      ForLoop< Destroy, 0, dim >::apply( mappings_ );
      integral_constant< int, 0 > codim0Variable;
      if( mappings_[ codim0Variable ].size() )
        delete mappings_[ codim0Variable ][ 0 ];
    }

    template< class Topology > class CornerStorage;
    template< int codim > struct Create;
    template< int codim > struct Destroy;

    struct GeometryTraits
      : public GenericGeometry::DefaultGeometryTraits< ctype, dim, dim >
    {
      typedef GenericGeometry::DefaultGeometryTraits< ctype, dim, dim > Base;

      typedef typename Base::CoordTraits CoordTraits;

      template< class Topology >
      struct Mapping
      {
        typedef GenericGeometry::CornerMapping< CoordTraits, Topology, dim, CornerStorage< Topology >, true > type;
      };

      struct Caching
      {
        static const GenericGeometry::EvaluationType evaluateJacobianTransposed = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateJacobianInverseTransposed = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateIntegrationElement = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateNormal = GenericGeometry::PreCompute;
      };

    };

  public:
    using Base::type;

    /** \brief Collection of types depending on the codimension */
    template< int codim >
    struct Codim
    {
      //! type of mapping embedding a subentity into the reference element
      typedef GenericGeometry::HybridMapping< dim-codim, GeometryTraits > Mapping;
    };

    /** \brief position of the barycenter of entity (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. This method returns the coordinates of
     *  the center of gravity of E within the current reference element.
     *
     *  \param[in]  i   number of subentity E (0 <= i < size( c ))
     *  \param[in]  c   codimension of subentity E
     */
    const FieldVector< ctype, dim > &position( int i, int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return baryCenters_[ c ][ i ];
    }

    /** \brief check if a coordinate is in the reference element
     *
     *  This method returns true if the given local coordinate is within this
     *  reference element.
     *
     *  \param[in]  local  coordinates of the point
     */
    bool checkInside ( const FieldVector< ctype, dim > &local ) const
    {
      return GenericGeometry::checkInside( type().id(), dim, local, 1e-12 );
    }

    /** \brief check if a local coordinate is in the reference element of
     *         the i-th subentity E with codimension c of the current
     *         reference element.
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method return true, if the given local
     *  coordinate is within the reference element for the entity E.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  local  coordinates of the point with respect to the
     *                     reference element of E
     *  \param[in]  i      number of subentity E (0 <= i < size( c ))
     */
    template< int codim >
    bool checkInside ( const FieldVector< ctype, dim-codim > &local, int i ) const
    {
      return GenericGeometry::checkInside( type( i, codim ).id(), dim-codim, local, 1e-12 );
    }

    /** \brief map a local coordinate on subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method maps a point within the reference
     *  element of E into the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  local  coordinates of the point with respect to the reference
     *                     element of E
     *  \param[in]  i      number of subentity E (0 <= i < size( c ))
     *  \param[in]  c      codimension of subentity E
     *
     *  \note The runtime argument c is redundant and must equal codim.
     *
     *  \note This method is just an alias for
     *  \code
     *  mapping< codim >( i ).global( local );
     *  \endcode
     */
    template< int codim >
    FieldVector< ctype, dim >
    global ( const FieldVector< ctype, dim-codim > &local, int i, int c ) const
    {
      if( c != codim )
        DUNE_THROW( Exception, "Local Coordinate Type does not correspond to codimension c." );
      assert( c == codim );
      return mapping< codim >( i ).global( local );
    }

    /** \brief map a local coordinate on subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method maps a point within the reference
     *  element of E into the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  local  coordinates of the point with respect to the reference
     *                     element of E
     *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
     *
     *  \note This method is just an alias for
     *  \code
     *  mapping< codim >( i ).global( local );
     *  \endcode
     */
    template< int codim >
    FieldVector< ctype, dim >
    global ( const FieldVector< ctype, dim-codim > &local, int i ) const
    {
      return mapping< codim >( i ).global( local );
    }

    /** \brief obtain the embedding of subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method returns a
     *  \ref Dune::GenericGeometry::HybridMapping HybridMapping that maps
     *  the reference element of E into the current reference element.
     *
     *  This method can be used in a GenericGeometry to represent subentities
     *  of the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
     */
    template< int codim >
    typename Codim< codim >::Mapping &mapping( int i ) const
    {
      integral_constant< int, codim > codimVariable;
      return *(mappings_[ codimVariable ][ i ]);
    }

    /** \brief obtain the volume of the reference element */
    ctype volume () const
    {
      return volume_;
    }

    /** \brief obtain the volume outer normal of the reference element
     *
     *  The volume outer normal is the outer normal whose length coincides
     *  with the face's volume.
     *
     *  \param[in]  face  index of the face, whose normal is desired
     */
    const FieldVector< ctype, dim > &volumeOuterNormal ( int face ) const
    {
      assert( (face >= 0) && (face < int( volumeNormals_.size())) );
      return volumeNormals_[ face ];
    }

    /** \brief initialize the reference element
     *
     *  \tparam  Topology  topology of the desired reference element
     *
     *  \note The dimension of the topology must match dim.
     */
    template< class Topology >
    void initializeTopology ()
    {
      dune_static_assert( (Topology::dimension == dim),
                          "Cannot initialize reference element for different dimension." );
      typedef GenericGeometry::ReferenceDomain< Topology > ReferenceDomain;

      const unsigned int topologyId = Topology::id;
      Base::template initializeTopology( topologyId );

      // compute corners
      const unsigned int numVertices = Base::size( dim );
      baryCenters_[ dim ].resize( numVertices );
      GenericGeometry::referenceCorners( topologyId, dim, &(baryCenters_[ dim ][ 0 ]) );

      // compute barycenters
      for( int codim = 0; codim < dim; ++codim )
      {
        const unsigned int size = Base::size( codim );
        baryCenters_[ codim ].resize( size );
        for( unsigned int i = 0; i < size; ++i )
        {
          baryCenters_[ codim ][ i ] = FieldVector< ctype, dim >( ctype( 0 ) );
          const unsigned int numCorners = Base::size( i, codim, dim );
          for( unsigned int j = 0; j < numCorners; ++j )
            baryCenters_[ codim ][ i ] += baryCenters_[ dim ][ Base::subEntity( i, codim, j, dim ) ];
          baryCenters_[ codim ][ i ] *= ctype( 1 ) / ctype( numCorners );
        }
      }

      // compute reference element volume
      volume_ = GenericGeometry::template referenceVolume< ctype >( topologyId, dim );

      // compute normals
      volumeNormals_.resize( ReferenceDomain::numNormals );
      for( unsigned int i = 0; i < ReferenceDomain::numNormals; ++i )
        ReferenceDomain::integrationOuterNormal( i ,volumeNormals_[ i ] );

      // set up mappings
      typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits > VirtualMapping;

      integral_constant< int, 0 > codim0Variable;
      mappings_[ codim0Variable ].resize( 1 );
      mappings_[ codim0Variable ][ 0 ]  = new VirtualMapping( codim0Variable );

      Dune::ForLoop< Create, 0, dim >::apply( static_cast< Base & >( *this ), mappings_ );
    }

  private:
    /** \brief Stores all subentities of a given codimension */
    template< int codim >
    struct MappingArray
      : public std::vector< typename Codim< codim >::Mapping * >
    {};

    /** \brief Type to store all subentities of all codimensions */
    typedef GenericGeometry::CodimTable< MappingArray, dim > MappingsTable;

    /** \brief The reference element volume */
    ctype volume_;

    std::vector< FieldVector< ctype, dim > > baryCenters_[ dim+1 ];
    std::vector< FieldVector< ctype, dim > > volumeNormals_;

    /** \brief Stores all subentities of all codimensions */
    MappingsTable mappings_;
  };



  // ReferenceElement::CornerStorage
  // -------------------------------

  template< class ctype, int dim >
  template< class Topology >
  class ReferenceElement< ctype, dim >::CornerStorage
  {
    typedef GenericGeometry::ReferenceDomain< Topology > RefDomain;

  public:
    static const unsigned int size = Topology::numCorners;

    template< class SubTopology >
    struct SubStorage
    {
      typedef CornerStorage< SubTopology > type;
    };

    explicit CornerStorage ( const integral_constant< int, 0 > & )
    {
      for( unsigned int i = 0; i < size; ++i )
        RefDomain::corner( i, coords_[ i ] );
    }

    template< class Mapping, unsigned int codim >
    explicit
    CornerStorage ( const GenericGeometry::SubMappingCoords< Mapping, codim > &coords )
    {
      for( unsigned int i = 0; i < size; ++i )
        coords_[ i ] = coords[ i ];
    }

    const FieldVector< ctype, dim > &operator[] ( unsigned int i ) const
    {
      return coords_[ i ];
    }

  private:
    FieldVector< ctype, dim > coords_[ size ];
  };



  // ReferenceElement::Create
  // ------------------------

  template< class ctype, int dim >
  template< int codim >
  struct ReferenceElement< ctype, dim >::Create
  {
    static void apply ( const ReferenceElement< void, dim > &refElement, MappingsTable &mappings )
    {
      if( codim > 0 )
      {
        integral_constant< int, 0 > codim0Variable;
        const typename Codim< 0 >::Mapping &refMapping = *(mappings[ codim0Variable ][ 0 ]);

        typedef typename GenericGeometry::MappingProvider< typename Codim< 0 >::Mapping, codim > MappingProvider;
        integral_constant< int, codim > codimVariable;

        const unsigned int size = refElement.size( codim );
        mappings[ codimVariable ].resize( size );
        for( unsigned int i = 0; i < size; ++i )
        {
          char *storage = new char[ MappingProvider::maxMappingSize ];
          mappings[ codimVariable ][ i ] = refMapping.template trace< codim >( i, storage );
        }
      }
    }
  };



  // ReferenceElement::Destroy
  // -------------------------

  template< class ctype, int dim >
  template< int codim >
  struct ReferenceElement< ctype, dim >::Destroy
  {
    static void apply ( MappingsTable &mappings )
    {
      if( codim > 0 )
      {
        integral_constant< int, codim > codimVariable;
        for( size_t i = 0; i < mappings[ codimVariable ].size(); ++i )
        {
          typedef typename Codim< codim >::Mapping Mapping;
          mappings[ codimVariable ][ i ]->~Mapping();
          char *storage = (char *)mappings[ codimVariable ][ i ];
          delete[]( storage );
        }
      }
    }
  };



  // ReferenceElementContainer
  // -------------------------

  template< class ctype, int dim >
  class ReferenceElementContainer
  {
    static const unsigned int numTopologies = (1u << dim);

  public:
    typedef ReferenceElement< ctype, dim > value_type;
    typedef const value_type *const_iterator;

    ReferenceElementContainer ()
    {
      ForLoop< Builder, 0, numTopologies-1 >::apply( values_ );
    }

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == dim );
      return values_[ type.id() ];
    }

    const value_type &simplex () const
    {
      return values_[ GenericGeometry::SimplexTopology< dim >::type::id ];
    }

    const value_type &cube () const
    {
      return values_[ GenericGeometry::CubeTopology< dim >::type::id ];
    }

    const value_type &pyramid () const
    {
      return values_[ GenericGeometry::PyramidTopology< dim >::type::id ];
    }

    const value_type &prism () const
    {
      return values_[ GenericGeometry::PrismTopology< dim >::type::id ];
    }

    const_iterator begin () const { return values_; }
    const_iterator end () const { return values_ + numTopologies; }

  private:
    template< int topologyId >
    struct Builder
    {
      static void apply ( value_type (&values)[ numTopologies ] )
      {
        typedef typename GenericGeometry::Topology< topologyId, dim >::type Topology;
        values[ topologyId ].template initializeTopology< Topology >();
      }
    };

    value_type values_[ numTopologies ];
  };



  // ReferenceElements
  // -----------------

  /** \brief Class providing access to the singletons of the
   *  reference elements. Special methods are available for
   *  simplex and cube elements of any dimension.
   *  The method general can be used to obtain the reference element
   *  for a given geometry type.
   *
   *  \ingroup GeometryReferenceElements
   */
  template< class ctype, int dim >
  struct ReferenceElements
  {
    typedef typename ReferenceElementContainer< ctype, dim >::const_iterator Iterator;

    //! get general reference elements
    static const ReferenceElement< ctype, dim > &
    general ( const GeometryType &type )
    {
      return container() ( type );
    }

    //! get simplex reference elements
    static const ReferenceElement< ctype, dim > &simplex ()
    {
      return container().simplex();
    }

    //! get hypercube reference elements
    static const ReferenceElement< ctype, dim > &cube ()
    {
      return container().cube();
    }

    static Iterator begin () { return container().begin(); }
    static Iterator end () { return container().end(); }

  private:
    static const ReferenceElementContainer< ctype, dim > &container ()
    {
      static ReferenceElementContainer< ctype, dim > container;
      return container;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
