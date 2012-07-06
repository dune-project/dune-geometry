// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_CACHED_MAPPING_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_CACHED_MAPPING_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/genericgeometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>
#include <dune/geometry/genericgeometry/mapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // Internal Forward Declarations
    // -----------------------------

    template< unsigned int, class >
    class CachedJacobianTransposed;

    template< unsigned int, class >
    class CachedJacobianInverseTransposed;



    // CachedStorage
    // -------------

    template< unsigned int dim, class GeometryTraits >
    class CachedStorage
    {
      friend class CachedJacobianTransposed< dim, GeometryTraits >;

    public:
      static const unsigned int dimension = dim;
      static const unsigned int dimWorld = GeometryTraits::dimWorld;

      typedef MappingTraits< typename GeometryTraits::CoordTraits, dimension, dimWorld > Traits;

      typedef typename GeometryTraits::Caching Caching;

      CachedStorage ()
        : affine( false ),
          jacobianTransposedComputed( false ),
          jacobianInverseTransposedComputed( false ),
          integrationElementComputed( false )
      {}

      typename Traits::JacobianTransposedType jacobianTransposed;
      typename Traits::JacobianType jacobianInverseTransposed;
      typename Traits::FieldType integrationElement;

      bool affine : 1;

      bool jacobianTransposedComputed : 1;        // = affine, if jacobian transposed was computed
      bool jacobianInverseTransposedComputed : 1; // = affine, if jacobian inverse transposed was computed
      bool integrationElementComputed : 1;        // = affine, if integration element was computed
    };



    // CachedJacobianTranposed
    // -----------------------

    template< unsigned int dim, class GeometryTraits >
    class CachedJacobianTransposed
    {
      friend class CachedJacobianInverseTransposed< dim, GeometryTraits >;

      typedef CachedStorage< dim, GeometryTraits > Storage;
      typedef typename Storage::Traits Traits;

      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      typedef typename Traits::FieldType ctype;

      static const int rows = Traits::dimension;
      static const int cols = Traits::dimWorld;

      typedef typename Traits::JacobianTransposedType FieldMatrix;

      operator bool () const
      {
        return storage().jacobianTransposedComputed;
      }

      operator const FieldMatrix & () const
      {
        return storage().jacobianTransposed;
      }

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mv( x, y );
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mtv( x, y );
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umv( x, y );
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umtv( x, y );
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmv( x, y );
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmtv( x, y );
      }

      ctype det () const
      {
        if( !storage().integrationElementComputed )
        {
          storage().integrationElement = MatrixHelper::template sqrtDetAAT< rows, cols >( storage().jacobianTransposed );
          storage().integrationElementComputed = storage().affine;
        }
        return storage().integrationElement;
      }

    private:
      Storage &storage () const { return storage_; }

      mutable Storage storage_;
    };



    // CachedJacobianInverseTransposed
    // -------------------------------

    template< unsigned int dim, class GeometryTraits >
    class CachedJacobianInverseTransposed
    {
      template< class, class > friend class CachedMapping;

      typedef CachedJacobianTransposed< dim, GeometryTraits > JacobianTransposed;
      typedef typename JacobianTransposed::Storage Storage;
      typedef typename JacobianTransposed::Traits Traits;

      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      typedef typename Traits::FieldType ctype;

      static const int rows = Traits::dimWorld;
      static const int cols = Traits::dimension;

      typedef typename Traits::JacobianType FieldMatrix;

      operator bool () const
      {
        return storage().jacobianInverseTransposedComputed;
      }

      operator const FieldMatrix & () const
      {
        return storage().jacobianInverseTransposed;
      }

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mv( x, y );
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mtv( x, y );
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umv( x, y );
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).umtv( x, y );
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmv( x, y );
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        static_cast< const FieldMatrix & >( *this ).mmtv( x, y );
      }

      ctype det () const
      {
        // integrationElement is always computed with jacobianInverseTransposed
        return ctype( 1 ) / storage().integrationElement;
      }

    private:
      JacobianTransposed &jacobianTransposed () { return jacobianTransposed_; }
      const JacobianTransposed &jacobianTransposed () const { return jacobianTransposed_; }

      Storage &storage () const { return jacobianTransposed().storage(); }

      JacobianTransposed jacobianTransposed_;
    };



    // CachedMapping
    // -------------

    template< class Topology, class GeometryTraits >
    class CachedMapping
    {
      typedef CachedMapping< Topology, GeometryTraits > This;

      typedef typename GeometryTraits::template Mapping< Topology >::type MappingImpl;

    public:
      typedef MappingTraits< typename GeometryTraits::CoordTraits, Topology::dimension, GeometryTraits::dimWorld > Traits;

      typedef GenericGeometry::Mapping< typename GeometryTraits::CoordTraits, Topology, GeometryTraits::dimWorld, MappingImpl > Mapping;

      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      typedef CachedStorage< dimension, GeometryTraits > Storage;
      typedef CachedJacobianTransposed< dimension, GeometryTraits > JacobianTransposed;
      typedef CachedJacobianInverseTransposed< dimension, GeometryTraits > JacobianInverseTransposed;

      typedef GenericGeometry::ReferenceElement< Topology, FieldType > ReferenceElement;

      // can we safely assume that this mapping is always affine?
      static const bool alwaysAffine = Mapping::alwaysAffine;

      typedef typename GeometryTraits::Caching Caching;

    private:
      typedef typename Traits::MatrixHelper MatrixHelper;

    public:
      template< class CoordVector >
      explicit CachedMapping ( const CoordVector &coords )
        : mapping_( coords )
      {
        if( alwaysAffine )
          storage().affine = true;
        else
          computeJacobianTransposed( baryCenter() );
        preCompute();
      }

      template< class CoordVector >
      explicit CachedMapping ( const std::pair< const CoordVector &, bool > &coords )
        : mapping_( coords.first )
      {
        storage().affine = coords.second;
        preCompute();
      }

      bool affine () const { return (alwaysAffine || storage().affine); }
      Dune::GeometryType type () const { return Dune::GeometryType( Topology() ); }

      int numCorners () const { return ReferenceElement::numCorners; }
      GlobalCoordinate corner ( int i ) const { return mapping().corner( i ); }
      GlobalCoordinate center () const { return global( ReferenceElement::baryCenter() ); }

      static bool checkInside ( const LocalCoordinate &x ) { return ReferenceElement::checkInside( x ); }

      GlobalCoordinate global ( const LocalCoordinate &x ) const
      {
        GlobalCoordinate y;
        if( jacobianTransposed() )
        {
          y = corner( 0 );
          jacobianTransposed().umtv( x, y );
          //MatrixHelper::template ATx< dimension, dimWorld >( jacobianTransposed_, x, y );
          //y += corner( 0 );
        }
        else
          mapping().global( x, y );
        return y;
      }

      LocalCoordinate local ( const GlobalCoordinate &y ) const
      {
        LocalCoordinate x;
        if( jacobianInverseTransposed() )
        {
          GlobalCoordinate z = y - corner( 0 );
          jacobianInverseTransposed().mtv( z, x );
          // MatrixHelper::template ATx< dimWorld, dimension >( jacobianInverseTransposed(), z, x );
        }
        else if( affine() )
        {
          const JacobianTransposed &JT = jacobianTransposed( baryCenter() );
          GlobalCoordinate z = y - corner( 0 );
          MatrixHelper::template xTRightInvA< dimension, dimWorld >( JT, z, x );
        }
        else
          mapping().local( y, x );
        return x;
      }

      FieldType integrationElement ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluateI = Caching::evaluateIntegrationElement;
        const EvaluationType evaluateJ = Caching::evaluateJacobianInverseTransposed;
        if( ((evaluateI == PreCompute) || (evaluateJ == PreCompute)) && alwaysAffine )
          return storage().integrationElement;
        else
          return jacobianTransposed( x ).det();
      }

      FieldType volume () const
      {
        // do we need a quadrature of higher order, here?
        const FieldType refVolume = ReferenceElement::volume();
        return refVolume * integrationElement( baryCenter() );
      }

      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluate = Caching::evaluateJacobianTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianTransposed();

        if( !jacobianTransposed() )
          computeJacobianTransposed( x );
        return jacobianTransposed();
      }

      const JacobianInverseTransposed &
      jacobianInverseTransposed ( const LocalCoordinate &x ) const
      {
        const EvaluationType evaluate = Caching::evaluateJacobianInverseTransposed;
        if( (evaluate == PreCompute) && alwaysAffine )
          return jacobianInverseTransposed();

        if( !jacobianInverseTransposed() )
          computeJacobianInverseTransposed( x );
        return jacobianInverseTransposed();
      }

      const Mapping &mapping () const { return mapping_; }

    private:
      static const LocalCoordinate &baryCenter ()
      {
        return ReferenceElement::baryCenter();
      }

      Storage &storage () const
      {
        return jacobianInverseTransposed().storage();
      }

      const JacobianTransposed &jacobianTransposed () const
      {
        return jacobianInverseTransposed().jacobianTransposed();
      }

      const JacobianInverseTransposed &jacobianInverseTransposed () const
      {
        return jacobianInverseTransposed_;
      }

      void preCompute ()
      {
        assert( affine() == mapping().jacobianTransposed( baryCenter(), storage().jacobianTransposed ) );
        if( !affine() )
          return;

        if( (Caching::evaluateJacobianTransposed == PreCompute) && !jacobianTransposed() )
          computeJacobianTransposed( baryCenter() );

        if( Caching::evaluateJacobianInverseTransposed == PreCompute )
          computeJacobianInverseTransposed( baryCenter() );
        else if( Caching::evaluateIntegrationElement == PreCompute )
          jacobianTransposed().det();
      }

      void computeJacobianTransposed ( const LocalCoordinate &x ) const
      {
        storage().affine = mapping().jacobianTransposed( x, storage().jacobianTransposed );
        storage().jacobianTransposedComputed = affine();
      }

      void computeJacobianInverseTransposed ( const LocalCoordinate &x ) const
      {
        storage().integrationElement
          = MatrixHelper::template rightInvA< dimension, dimWorld >( jacobianTransposed( x ), storage().jacobianInverseTransposed );
        storage().integrationElementComputed = affine();
        storage().jacobianInverseTransposedComputed = affine();
      }

    private:
      Mapping mapping_;
      JacobianInverseTransposed jacobianInverseTransposed_;
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_CACHED_MAPPING_HH
