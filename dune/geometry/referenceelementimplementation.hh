// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_REFERENCEELEMENTIMPLEMENTATION_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTIMPLEMENTATION_HH

#include <cassert>

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include <array>
#include <bitset>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

namespace Dune
{

  namespace Geo
  {

#ifndef DOXYGEN

    // Internal Forward Declarations
    // -----------------------------

    namespace Impl
    {
      template< class ctype, int dim >
      class ReferenceElementContainer;
    }

    template< class ctype, int dim >
    struct ReferenceElements;



    namespace Impl
    {

      using Dune::Impl::isPrism;
      using Dune::Impl::isPyramid;
      using Dune::Impl::baseTopologyId;
      using Dune::Impl::prismConstruction;
      using Dune::Impl::pyramidConstruction;
      using Dune::Impl::numTopologies;

      /** \brief Compute the number of subentities of a given codimension */
      unsigned int size ( unsigned int topologyId, int dim, int codim );



      /** \brief Compute the topology id of a given subentity
       *
       * \param topologyId Topology id of entity
       * \param dim Dimension of entity
       * \param codim Codimension of the subentity that we are interested in
       * \param i Number of the subentity that we are interested in
       */
      unsigned int subTopologyId ( unsigned int topologyId, int dim, int codim, unsigned int i );



      // subTopologyNumbering
      // --------------------

      void subTopologyNumbering ( unsigned int topologyId, int dim, int codim, unsigned int i, int subcodim,
                                  unsigned int *beginOut, unsigned int *endOut );




      // checkInside
      // -----------

      template< class ct, int cdim >
      inline bool
      checkInside ( unsigned int topologyId, int dim, const FieldVector< ct, cdim > &x, ct tolerance, ct factor = ct( 1 ) )
      {
        assert( (dim >= 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( dim > 0 )
          {
            const ct baseFactor = (isPrism( topologyId, dim ) ? factor : factor - x[ dim-1 ]);
            if( (x[ dim-1 ] > -tolerance) && (factor - x[ dim-1 ] > -tolerance) )
              return checkInside< ct, cdim >( baseTopologyId( topologyId, dim ), dim-1, x, tolerance, baseFactor );
            else
              return false;
          }
        else
          return true;
      }



      // referenceCorners
      // ----------------

      template< class ct, int cdim >
      inline unsigned int
      referenceCorners ( unsigned int topologyId, int dim, FieldVector< ct, cdim > *corners )
      {
        assert( (dim >= 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( dim > 0 )
          {
            const unsigned int nBaseCorners
              = referenceCorners( baseTopologyId( topologyId, dim ), dim-1, corners );
            assert( nBaseCorners == size( baseTopologyId( topologyId, dim ), dim-1, dim-1 ) );
            if( isPrism( topologyId, dim ) )
              {
                std::copy( corners, corners + nBaseCorners, corners + nBaseCorners );
                for( unsigned int i = 0; i < nBaseCorners; ++i )
                  corners[ i+nBaseCorners ][ dim-1 ] = ct( 1 );
                return 2*nBaseCorners;
              }
            else
              {
                corners[ nBaseCorners ] = FieldVector< ct, cdim >( ct( 0 ) );
                corners[ nBaseCorners ][ dim-1 ] = ct( 1 );
                return nBaseCorners+1;
              }
          }
        else
          {
            *corners = FieldVector< ct, cdim >( ct( 0 ) );
            return 1;
          }
      }



      // referenceVolume
      // ---------------

      unsigned long referenceVolumeInverse ( unsigned int topologyId, int dim );

      template< class ct >
      inline ct referenceVolume ( unsigned int topologyId, int dim )
      {
        return ct( 1 ) / ct( referenceVolumeInverse( topologyId, dim ) );
      }



      // referenceOrigins
      // ----------------

      template< class ct, int cdim >
      inline unsigned int
      referenceOrigins ( unsigned int topologyId, int dim, int codim, FieldVector< ct, cdim > *origins )
      {
        assert( (dim >= 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );
        assert( (codim >= 0) && (codim <= dim) );

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? referenceOrigins( baseId, dim-1, codim, origins ) : 0);
                const unsigned int m = referenceOrigins( baseId, dim-1, codim-1, origins+n );
                for( unsigned int i = 0; i < m; ++i )
                  {
                    origins[ n+m+i ] = origins[ n+i ];
                    origins[ n+m+i ][ dim-1 ] = ct( 1 );
                  }
                return n+2*m;
              }
            else
              {
                const unsigned int m = referenceOrigins( baseId, dim-1, codim-1, origins );
                if( codim == dim )
                  {
                    origins[ m ] = FieldVector< ct, cdim >( ct( 0 ) );
                    origins[ m ][ dim-1 ] = ct( 1 );
                    return m+1;
                  }
                else
                  return m+referenceOrigins( baseId, dim-1, codim, origins+m );
              }
          }
        else
          {
            origins[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
            return 1;
          }
      }



      // referenceEmbeddings
      // -------------------

      template< class ct, int cdim, int mydim >
      inline unsigned int
      referenceEmbeddings ( unsigned int topologyId, int dim, int codim,
                            FieldVector< ct, cdim > *origins,
                            FieldMatrix< ct, mydim, cdim > *jacobianTransposeds )
      {
        assert( (0 <= codim) && (codim <= dim) && (dim <= cdim) );
        assert( (dim - codim <= mydim) && (mydim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? referenceEmbeddings( baseId, dim-1, codim, origins, jacobianTransposeds ) : 0);
                for( unsigned int i = 0; i < n; ++i )
                  jacobianTransposeds[ i ][ dim-codim-1 ][ dim-1 ] = ct( 1 );

                const unsigned int m = referenceEmbeddings( baseId, dim-1, codim-1, origins+n, jacobianTransposeds+n );
                std::copy( origins+n, origins+n+m, origins+n+m );
                std::copy( jacobianTransposeds+n, jacobianTransposeds+n+m, jacobianTransposeds+n+m );
                for( unsigned int i = 0; i < m; ++i )
                  origins[ n+m+i ][ dim-1 ] = ct( 1 );

                return n+2*m;
              }
            else
              {
                const unsigned int m = referenceEmbeddings( baseId, dim-1, codim-1, origins, jacobianTransposeds );
                if( codim == dim )
                  {
                    origins[ m ] = FieldVector< ct, cdim >( ct( 0 ) );
                    origins[ m ][ dim-1 ] = ct( 1 );
                    jacobianTransposeds[ m ] = FieldMatrix< ct, mydim, cdim >( ct( 0 ) );
                    return m+1;
                  }
                else
                  {
                    const unsigned int n = referenceEmbeddings( baseId, dim-1, codim, origins+m, jacobianTransposeds+m );
                    for( unsigned int i = 0; i < n; ++i )
                      {
                        for( int k = 0; k < dim-1; ++k )
                          jacobianTransposeds[ m+i ][ dim-codim-1 ][ k ] = -origins[ m+i ][ k ];
                        jacobianTransposeds[ m+i ][ dim-codim-1 ][ dim-1 ] = ct( 1 );
                      }
                    return m+n;
                  }
              }
          }
        else
          {
            origins[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
            jacobianTransposeds[ 0 ] = FieldMatrix< ct, mydim, cdim >( ct( 0 ) );
            for( int k = 0; k < dim; ++k )
              jacobianTransposeds[ 0 ][ k ][ k ] = ct( 1 );
            return 1;
          }
      }



      // referenceIntegrationOuterNormals
      // --------------------------------

      template< class ct, int cdim >
      inline unsigned int
      referenceIntegrationOuterNormals ( unsigned int topologyId, int dim,
                                         const FieldVector< ct, cdim > *origins,
                                         FieldVector< ct, cdim > *normals )
      {
        assert( (dim > 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( dim > 1 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            if( isPrism( topologyId, dim ) )
              {
                const unsigned int numBaseFaces
                  = referenceIntegrationOuterNormals( baseId, dim-1, origins, normals );

                for( unsigned int i = 0; i < 2; ++i )
                  {
                    normals[ numBaseFaces+i ] = FieldVector< ct, cdim >( ct( 0 ) );
                    normals[ numBaseFaces+i ][ dim-1 ] = ct( 2*int( i )-1 );
                  }

                return numBaseFaces+2;
              }
            else
              {
                normals[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
                normals[ 0 ][ dim-1 ] = ct( -1 );

                const unsigned int numBaseFaces
                  = referenceIntegrationOuterNormals( baseId, dim-1, origins+1, normals+1 );
                for( unsigned int i = 1; i <= numBaseFaces; ++i )
                  normals[ i ][ dim-1 ] = normals[ i ]*origins[ i ];

                return numBaseFaces+1;
              }
          }
        else
          {
            for( unsigned int i = 0; i < 2; ++i )
              {
                normals[ i ] = FieldVector< ct, cdim >( ct( 0 ) );
                normals[ i ][ 0 ] = ct( 2*int( i )-1 );
              }

            return 2;
          }
      }

      template< class ct, int cdim >
      inline unsigned int
      referenceIntegrationOuterNormals ( unsigned int topologyId, int dim,
                                         FieldVector< ct, cdim > *normals )
      {
        assert( (dim > 0) && (dim <= cdim) );

        FieldVector< ct, cdim > *origins
          = new FieldVector< ct, cdim >[ size( topologyId, dim, 1 ) ];
        referenceOrigins( topologyId, dim, 1, origins );

        const unsigned int numFaces
          = referenceIntegrationOuterNormals( topologyId, dim, origins, normals );
        assert( numFaces == size( topologyId, dim, 1 ) );

        delete[] origins;

        return numFaces;
      }

    } // namespace Impl



    // ReferenceElement
    // ----------------

    /** \class ReferenceElementImplementation
     *  \ingroup GeometryReferenceElements
     *  \brief This class provides access to geometric and topological
     *  properties of a reference element.
     *
     *  This includes its type,
     *  the number of subentities, the volume, and a method for checking
     *  if a point is contained in the reference element.
     *  The embedding of each subentity into the reference element is also
     *  provided.
     *
     *  A singleton of this class for a given geometry type can be accessed
     *  through the ReferenceElements class.

     *  \tparam ctype  field type for coordinates
     *  \tparam dim    dimension of the reference element
     *
     */
    template< class ctype_, int dim >
    class ReferenceElementImplementation
    {

    public:

      //! The coordinate field type.
      using ctype = ctype_;

      //! The coordinate field type.
      using CoordinateField = ctype;

      //! The coordinate type.
      using Coordinate = Dune::FieldVector<ctype,dim>;

      //! The dimension of the reference element.
      static constexpr int dimension = dim;

      /** \brief Type used for volume */
      typedef ctype Volume;

    private:

      friend class Impl::ReferenceElementContainer< ctype, dim >;

      struct SubEntityInfo;

      template< int codim > struct CreateGeometries;

    public:
      /** \brief Collection of types depending on the codimension */
      template< int codim >
      struct Codim
      {
        //! type of geometry embedding a subentity into the reference element
        typedef AffineGeometry< ctype, dim-codim, dim > Geometry;
      };

      // ReferenceElement cannot be copied.
      ReferenceElementImplementation ( const ReferenceElementImplementation& ) = delete;

      // ReferenceElementImplementation cannot be copied.
      ReferenceElementImplementation& operator= ( const ReferenceElementImplementation& ) = delete;

      // ReferenceElementImplementation is default-constructible (required for storage in std::array)
      ReferenceElementImplementation () = default;

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
       *  a subentity of E. If cc<c this number is zero.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E (0 <= c <= dim)
       *  \param[in]  cc  codimension whose size is desired (0 <= cc <= dim)
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
       *  (cc-c) of E. Then, S is a also a subentity of codimension cc of the current
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

      /** \brief Obtain the range of numbers of subentities with codim cc of (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns a range of numbers of
       *  all subentities of E with codimension cc. Notice that the sub-subentity
       *  codimension as well as the numbers in the returned range are
       *  given with respect to the reference element itself and not with
       *  respect to E. For 0<=cc<c this will return an empty range.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  cc  codimension of subentity S (0 <= cc <= dim)
       *
       *  \returns An iterable range of numbers of the sub-subentities.
       */
      auto subEntities ( int i, int c, int cc ) const
      {
        assert( (i >= 0) && (i < size( c )) );
        return info_[ c ][ i ].numbers( cc );
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
      const Coordinate &position( int i, int c ) const
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
      bool checkInside ( const Coordinate &local ) const
      {
        const ctype tolerance = ctype( 64 ) * std::numeric_limits< ctype >::epsilon();
        return Impl::template checkInside< ctype, dim >( type().id(), dim, local, tolerance );
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
        return std::get< codim >( geometries_ )[ i ];
      }

      /** \brief obtain the volume of the reference element */
      Volume volume () const
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
      const Coordinate &integrationOuterNormal ( int face ) const
      {
        assert( (face >= 0) && (face < int( integrationNormals_.size() )) );
        return integrationNormals_[ face ];
      }

    private:
      void initialize ( unsigned int topologyId )
      {
        assert( topologyId < Impl::numTopologies( dim ) );

        // set up subentities
        for( int codim = 0; codim <= dim; ++codim )
          {
            const unsigned int size = Impl::size( topologyId, dim, codim );
            info_[ codim ].resize( size );
            for( unsigned int i = 0; i < size; ++i )
              info_[ codim ][ i ].initialize( topologyId, codim, i );
          }

        // compute corners
        const unsigned int numVertices = size( dim );
        baryCenters_[ dim ].resize( numVertices );
        Impl::referenceCorners( topologyId, dim, &(baryCenters_[ dim ][ 0 ]) );

        // compute barycenters
        for( int codim = 0; codim < dim; ++codim )
          {
            baryCenters_[ codim ].resize( size(codim) );
            for( int i = 0; i < size( codim ); ++i )
              {
                baryCenters_[ codim ][ i ] = Coordinate( ctype( 0 ) );
                const unsigned int numCorners = size( i, codim, dim );
                for( unsigned int j = 0; j < numCorners; ++j )
                  baryCenters_[ codim ][ i ] += baryCenters_[ dim ][ subEntity( i, codim, j, dim ) ];
                baryCenters_[ codim ][ i ] *= ctype( 1 ) / ctype( numCorners );
              }
          }

        // compute reference element volume
        volume_ = Impl::template referenceVolume< ctype >( topologyId, dim );

        // compute integration outer normals
        if( dim > 0 )
          {
            integrationNormals_.resize( size( 1 ) );
            Impl::referenceIntegrationOuterNormals( topologyId, dim, &(integrationNormals_[ 0 ]) );
          }

        // set up geometries
        Hybrid::forEach( std::make_index_sequence< dim+1 >{}, [ & ]( auto i ){ CreateGeometries< i >::apply( *this, geometries_ ); } );
      }

      template< int... codim >
      static std::tuple< std::vector< typename Codim< codim >::Geometry >... >
      makeGeometryTable ( std::integer_sequence< int, codim... > );

      /** \brief Type to store all subentities of all codimensions */
      typedef decltype( makeGeometryTable( std::make_integer_sequence< int, dim+1 >() ) ) GeometryTable;

      /** \brief The reference element volume */
      ctype volume_;

      std::vector< Coordinate > baryCenters_[ dim+1 ];
      std::vector< Coordinate > integrationNormals_;

      /** \brief Stores all subentities of all codimensions */
      GeometryTable geometries_;

      std::vector< SubEntityInfo > info_[ dim+1 ];
    };

    /** \brief topological information about the subentities of a reference element */
    template< class ctype, int dim >
    struct ReferenceElementImplementation< ctype, dim >::SubEntityInfo
    {
      // Compute upper bound for the number of subsentities.
      // If someone knows an explicit formal feel free to
      // implement it here.
      static constexpr std::size_t maxSubEntityCount()
      {
        std::size_t maxCount=0;
        for(std::size_t codim=0; codim<=dim; ++codim)
          maxCount = std::max(maxCount, binomial(std::size_t(dim),codim)*(1 << codim));
        return maxCount;
      }

      using SubEntityFlags = std::bitset<maxSubEntityCount()>;

      class SubEntityRange
        : public Dune::IteratorRange<const unsigned int*>
      {
        using Base = typename Dune::IteratorRange<const unsigned int*>;

      public:

        using iterator = Base::iterator;
        using const_iterator = Base::const_iterator;

        SubEntityRange(const iterator& begin, const iterator& end, const SubEntityFlags& contains) :
          Base(begin, end),
          containsPtr_(&contains),
          size_(end-begin)
        {}

        SubEntityRange() :
          Base(),
          containsPtr_(nullptr),
          size_(0)
        {}

        std::size_t size() const
        {
          return size_;
        }

        bool contains(std::size_t i) const
        {
          return (*containsPtr_)[i];
        }

      private:
        const SubEntityFlags* containsPtr_;
        std::size_t size_;
        std::size_t offset_;
      };

      using NumberRange = typename Dune::IteratorRange<const unsigned int*>;

      SubEntityInfo ()
        : numbering_( nullptr )
      {
        std::fill( offset_.begin(), offset_.end(), 0 );
      }

      SubEntityInfo ( const SubEntityInfo &other )
        : offset_( other.offset_ ),
          type_( other.type_ ),
          containsSubentity_( other.containsSubentity_ )
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

        containsSubentity_ = other.containsSubentity_;

        return *this;
      }

      int size ( int cc ) const
      {
        assert( (cc >= 0) && (cc <= dim) );
        return (offset_[ cc+1 ] - offset_[ cc ]);
      }

      int number ( int ii, int cc ) const
      {
        assert( (ii >= 0) && (ii < size( cc )) );
        return numbering_[ offset_[ cc ] + ii ];
      }

      auto numbers ( int cc ) const
      {
        return SubEntityRange( numbering_ + offset_[ cc ], numbering_ + offset_[ cc+1 ], containsSubentity_[cc]);
      }

      const GeometryType &type () const { return type_; }

      void initialize ( unsigned int topologyId, int codim, unsigned int i )
      {
        const unsigned int subId = Impl::subTopologyId( topologyId, dim, codim, i );
        type_ = GeometryType( subId, dim-codim );

        // compute offsets
        for( int cc = 0; cc <= codim; ++cc )
          offset_[ cc ] = 0;
        for( int cc = codim; cc <= dim; ++cc )
          offset_[ cc+1 ] = offset_[ cc ] + Impl::size( subId, dim-codim, cc-codim );

        // compute subnumbering
        deallocate( numbering_ );
        numbering_ = allocate();
        for( int cc = codim; cc <= dim; ++cc )
          Impl::subTopologyNumbering( topologyId, dim, codim, i, cc-codim, numbering_+offset_[ cc ], numbering_+offset_[ cc+1 ] );

        // initialize containsSubentity lookup-table
        for(std::size_t cc=0; cc<= dim; ++cc)
        {
          containsSubentity_[cc].reset();
          for(std::size_t idx=0; idx<std::size_t(size(cc)); ++idx)
            containsSubentity_[cc][number(idx,cc)] = true;
        }
      }

    protected:
      int codim () const { return dim - type().dim(); }

      unsigned int *allocate () { return (capacity() != 0 ? new unsigned int[ capacity() ] : nullptr); }
      void deallocate ( unsigned int *ptr ) { delete[] ptr; }
      unsigned int capacity () const { return offset_[ dim+1 ]; }

    private:
      unsigned int *numbering_;
      std::array< unsigned int, dim+2 > offset_;
      GeometryType type_;
      std::array< SubEntityFlags, dim+1> containsSubentity_;
    };


    template< class ctype, int dim >
    template< int codim >
    struct ReferenceElementImplementation< ctype, dim >::CreateGeometries
    {
      template< int cc >
      static typename ReferenceElements< ctype, dim-cc >::ReferenceElement
      subRefElement( const ReferenceElementImplementation< ctype, dim > &refElement, int i, std::integral_constant< int, cc > )
      {
        return ReferenceElements< ctype, dim-cc >::general( refElement.type( i, cc ) );
      }

      static typename ReferenceElements< ctype, dim >::ReferenceElement
      subRefElement(const ReferenceElementImplementation< ctype, dim > &refElement,
                    [[maybe_unused]] int i, std::integral_constant<int, 0>)
      {
        return refElement;
      }

      static void apply ( const ReferenceElementImplementation< ctype, dim > &refElement, GeometryTable &geometries )
      {
        const int size = refElement.size( codim );
        std::vector< FieldVector< ctype, dim > > origins( size );
        std::vector< FieldMatrix< ctype, dim - codim, dim > > jacobianTransposeds( size );
        Impl::referenceEmbeddings( refElement.type().id(), dim, codim, &(origins[ 0 ]), &(jacobianTransposeds[ 0 ]) );

        std::get< codim >( geometries ).reserve( size );
        for( int i = 0; i < size; ++i )
          {
            typename Codim< codim >::Geometry geometry( subRefElement( refElement, i, std::integral_constant< int, codim >() ), origins[ i ], jacobianTransposeds[ i ] );
            std::get< codim >( geometries ).push_back( geometry );
          }
      }
    };

#endif // DOXYGEN

  } // namespace Geo

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTIMPLEMENTATION_HH
