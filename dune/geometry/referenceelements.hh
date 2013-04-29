// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTS_HH

#include <algorithm>
#include <limits>

#include <dune/common/deprecated.hh>
#include <dune/common/array.hh>
#include <dune/common/forloop.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/genericgeometry/codimtable.hh>
#include <dune/geometry/genericgeometry/subtopologies.hh>
#include <dune/geometry/genericgeometry/referencedomain.hh>

// for backward compatibility include header
// with deprecated classes. Remove after Dune 2.3
#include "genericreferenceelements.hh"

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElementContainer;

  template< class ctype, int dim >
  struct ReferenceElements;



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
  {
    typedef ReferenceElement< ctype, dim > This;

    friend class ReferenceElementContainer< ctype, dim >;

    struct SubEntityInfo;

    // make copy constructor private
    ReferenceElement ( const This & );

    ReferenceElement () {}

    template< int codim > struct CreateGeometries;

  public:
    /** \brief Collection of types depending on the codimension */
    template< int codim >
    struct Codim
    {
      //! type of geometry embedding a subentity into the reference element
      typedef AffineGeometry< ctype, dim-codim, dim > Geometry;
      typedef Geometry Mapping DUNE_DEPRECATED_MSG ( "Use Geometry instead." );
    };

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
      const ctype tolerance = ctype( 64 ) * std::numeric_limits< ctype >::epsilon();
      return GenericGeometry::template checkInside< ctype, dim >( type().id(), dim, local, tolerance );
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
      const ctype tolerance = ctype( 64 ) * std::numeric_limits< ctype >::epsilon();
      return GenericGeometry::template checkInside< ctype, dim >( type( i, codim ).id(), dim-codim, local, tolerance );
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
     */
    template< int codim >
    FieldVector< ctype, dim >
    DUNE_DEPRECATED_MSG( "Use geometry< codim >( i ).global( local ) instead." )
    global( const FieldVector< ctype, dim-codim > &local, int i, int c ) const
    {
      if( c != codim )
        DUNE_THROW( Exception, "Local Coordinate Type does not correspond to codimension c." );
      assert( c == codim );
      return geometry< codim >( i ).global( local );
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
     */
    template< int codim >
    FieldVector< ctype, dim >
    DUNE_DEPRECATED_MSG( "Use geometry< codim >( i ).global( local ) instead." )
    global( const FieldVector< ctype, dim-codim > &local, int i ) const
    {
      return geometry< codim >( i ).global( local );
    }

    /** \brief obtain the embedding of subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method returns a \ref Dune::AffineGeometry
     *  that maps the reference element of E into the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
     */
    template< int codim >
    typename Codim< codim >::Geometry geometry ( int i ) const
    {
      integral_constant< int, codim > codimVariable;
      return geometries_[ codimVariable ][ i ];
    }

    /** \brief obtain the embedding of subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method returns a \ref Dune::AffineGeometry
     *  that maps the reference element of E into the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
     */
    template< int codim >
    DUNE_DEPRECATED_MSG( "Use geometry(i) instead." )
    const typename Codim< codim >::Mapping &mapping( int i ) const
    {
      integral_constant< int, codim > codimVariable;
      return geometries_[ codimVariable ][ i ];
    }

    /** \brief obtain the volume of the reference element */
    ctype volume () const
    {
      return volume_;
    }

    /** \brief obtain the integration outer normal of the reference element
     *
     *  The integration outer normal is the outer normal whose length coincides
     *  with the face's integration element.
     *
     *  \param[in]  face  index of the face, whose normal is desired
     */
    const FieldVector< ctype, dim > &integrationOuterNormal ( int face ) const
    {
      assert( (face >= 0) && (face < int( integrationNormals_.size() )) );
      return integrationNormals_[ face ];
    }

    const FieldVector< ctype, dim > &volumeOuterNormal ( int face ) const
    DUNE_DEPRECATED_MSG( "This method has always returned the integration outer normal; use integrationOuterNormal instead." )
    {
      return integrationOuterNormal( face );
    }

    /** \brief initialize the reference element
     *
     *  \param[in]  topologyId  topology id for the desired reference element
     */
    void initializeTopology ( unsigned int topologyId )
    DUNE_DEPRECATED_MSG( "initializeTopology has never been an official interface method." )
    {
      initialize( topologyId );
    }

  private:
    void initialize ( unsigned int topologyId )
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

      // compute corners
      const unsigned int numVertices = size( dim );
      baryCenters_[ dim ].resize( numVertices );
      GenericGeometry::referenceCorners( topologyId, dim, &(baryCenters_[ dim ][ 0 ]) );

      // compute barycenters
      for( int codim = 0; codim < dim; ++codim )
      {
        baryCenters_[ codim ].resize( size(codim) );
        for( int i = 0; i < size( codim ); ++i )
        {
          baryCenters_[ codim ][ i ] = FieldVector< ctype, dim >( ctype( 0 ) );
          const unsigned int numCorners = size( i, codim, dim );
          for( unsigned int j = 0; j < numCorners; ++j )
            baryCenters_[ codim ][ i ] += baryCenters_[ dim ][ subEntity( i, codim, j, dim ) ];
          baryCenters_[ codim ][ i ] *= ctype( 1 ) / ctype( numCorners );
        }
      }

      // compute reference element volume
      volume_ = GenericGeometry::template referenceVolume< ctype >( topologyId, dim );

      // compute integration outer normals
      if( dim > 0 )
      {
        integrationNormals_.resize( size( 1 ) );
        GenericGeometry::referenceIntegrationOuterNormals( topologyId, dim, &(integrationNormals_[ 0 ]) );
      }

      // set up geometries
      Dune::ForLoop< CreateGeometries, 0, dim >::apply( *this, geometries_ );
    }

    /** \brief Stores all subentities of a given codimension */
    template< int codim >
    struct GeometryArray
    : public std::vector< typename Codim< codim >::Geometry >
    {};

    /** \brief Type to store all subentities of all codimensions */
    typedef GenericGeometry::CodimTable< GeometryArray, dim > GeometryTable;

    /** \brief The reference element volume */
    ctype volume_;

    std::vector< FieldVector< ctype, dim > > baryCenters_[ dim+1 ];
    std::vector< FieldVector< ctype, dim > > integrationNormals_;

    /** \brief Stores all subentities of all codimensions */
    GeometryTable geometries_;

    std::vector< SubEntityInfo > info_[ dim+1 ];
  };

  /** \brief topological information about the subentities of a reference element */
  template< class ctype, int dim >
  struct ReferenceElement< ctype, dim >::SubEntityInfo
  {
    SubEntityInfo ()
      : numbering_( nullptr )
    {
      std::fill( offset_.begin(), offset_.end(), 0 );
    }

    SubEntityInfo ( const SubEntityInfo &other )
      : offset_( other.offset_ ),
        type_( other.type_ )
    {
      numbering_ = allocate();
      std::copy( other.numbering_, other.numbering_ + capacity(), numbering_ );
    }

    ~SubEntityInfo () { deallocate( numbering_ ); }

    const SubEntityInfo &operator= ( const SubEntityInfo &other )
    {
      type_ = other.type_;
      offset_ = other.offset_;

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
    array< unsigned int, dim+2 > offset_;
    GeometryType type_;
  };


  template< class ctype, int dim >
  template< int codim >
  struct ReferenceElement< ctype, dim >::CreateGeometries
  {
    template< int cc >
    static const ReferenceElement< ctype, dim-cc > &
    subRefElement( const ReferenceElement< ctype, dim > &refElement, int i, integral_constant< int, cc > )
    {
      return ReferenceElements< ctype, dim-cc >::general( refElement.type( i, cc ) );
    }

    static const ReferenceElement< ctype, dim > &
    subRefElement( const ReferenceElement< ctype, dim > &refElement, int i, integral_constant< int, 0 > )
    {
      return refElement;
    }

    static void apply ( const ReferenceElement< ctype, dim > &refElement, GeometryTable &geometries )
    {
      const int size = refElement.size( codim );
      std::vector< FieldVector< ctype, dim > > origins( size );
      std::vector< FieldMatrix< ctype, dim - codim, dim > > jacobianTransposeds( size );
      GenericGeometry::referenceEmbeddings( refElement.type().id(), dim, codim, &(origins[ 0 ]), &(jacobianTransposeds[ 0 ]) );

      integral_constant< int, codim > codimVariable;
      geometries[ codimVariable ].reserve( size );
      for( int i = 0; i < size; ++i )
      {
        typename Codim< codim >::Geometry geometry( subRefElement( refElement, i, codimVariable ), origins[ i ], jacobianTransposeds[ i ] );
        geometries[ codimVariable ].push_back( geometry );
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
      for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        values_[ topologyId ].initialize( topologyId );
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
    value_type values_[ numTopologies ];
  };



  // ReferenceElements
  // ------------------------

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
    DUNE_EXPORT static const ReferenceElementContainer< ctype, dim > &container ()
    {
      static ReferenceElementContainer< ctype, dim > container;
      return container;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
