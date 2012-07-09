// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_HYBRIDMAPPING_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_HYBRIDMAPPING_HH

#include <cstddef>

#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/cachedmapping.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/traceprovider.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // Internal Forward Declarations
    // -----------------------------

    template< unsigned int dim, class GeometryTraits >
    class HybridMapping;

    template< class Topology, class GeometryTraits >
    class VirtualMapping;



    // HybridMappingBase
    // -----------------

    /** \cond */
    template< unsigned int dim, class GeometryTraits, unsigned int codim = dim >
    class HybridMappingBase;

    template< unsigned int dim, class GeometryTraits, unsigned int codim >
    class HybridMappingBase
      : public virtual HybridMappingBase< dim, GeometryTraits, codim-1 >
    {
      typedef HybridMapping< dim, GeometryTraits > Mapping;

    public:
      virtual ~HybridMappingBase() {}

    protected:
      using HybridMappingBase< dim, GeometryTraits, codim-1 >::trace;

      virtual HybridMapping< dim - codim, GeometryTraits > *
      trace ( integral_constant< int, codim >, unsigned int i, char *mappingStorage ) const = 0;
    };

    template< unsigned int dim, class GeometryTraits >
    class HybridMappingBase< dim, GeometryTraits, 0 >
    {
      typedef HybridMapping< dim, GeometryTraits > Mapping;

    public:
      virtual ~HybridMappingBase() {}

    protected:
      virtual HybridMapping< dim, GeometryTraits > *
      trace ( integral_constant< int, 0 >, unsigned int i, char *mappingStorage ) const = 0;
    };
    /** \endcond */



    // HybridMapping
    // -------------

    /** \class   HybridMapping
     *  \ingroup GenericGeometry
     *  \brief   abstract base class for generic mapping
     *
     *  This is the user-visible class of the generic geometries if the
     *  topology type for each codimension is not unique.
     */
    template< unsigned int dim, class GeometryTraits >
    class HybridMapping
    /** \cond */
      : public virtual HybridMappingBase< dim, GeometryTraits >
        /** \endcond */
    {
      typedef HybridMapping< dim, GeometryTraits > This;

    protected:
      typedef MappingTraits< typename GeometryTraits::CoordTraits, dim, GeometryTraits::dimWorld > Traits;

    public:
      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      typedef CachedJacobianTransposed< dimension, GeometryTraits > JacobianTransposed;
      typedef CachedJacobianInverseTransposed< dimension, GeometryTraits > JacobianInverseTransposed;

      template< int codim >
      struct Codim
      {
        typedef HybridMapping< dimension - codim, GeometryTraits > Trace;
      };

      typedef typename GeometryTraits::Caching Caching;
      typedef typename GeometryTraits::UserData UserData;

      virtual ~HybridMapping ()
      {}

      /** \brief is this mapping affine? */
      virtual bool affine () const = 0;
      /** \brief obtain the name of the reference element */
      virtual Dune::GeometryType type () const = 0;

      /** \brief obtain number of corners of the corresponding reference element */
      virtual int numCorners () const = 0;
      /** \brief obtain coordinates of the i-th corner */
      virtual GlobalCoordinate corner ( int i ) const = 0;
      /** \brief obtain the centroid of the mapping's image */
      virtual GlobalCoordinate center () const = 0;

      /** \brief evaluate the mapping
       *
       *  \param[in]  x  local coordinate to map
       *
       *  \returns corresponding global coordinate
       */
      virtual GlobalCoordinate global ( const LocalCoordinate &x ) const = 0;
      /** \brief evaluate the inverse mapping
       *
       *  \param[in]  y  global coordinate to map
       *
       *  \return corresponding local coordinate
       *
       *  \note The returned local coordinate y minimizes
       *  \code
       *  (global( x ) - y).two_norm()
       *  \endcode
       */
      virtual LocalCoordinate local ( const GlobalCoordinate &y ) const = 0;

      /** \brief check whether a point lies within the reference element
       *
       *  \param[in]  x  local coordinate of point to check
       *
       *  \note Historically, this method was part of the geometry interface.
       *        It is still required for the GenericReferenceElement.
       */
      virtual bool checkInside ( const LocalCoordinate &x ) const = 0;

      /** \brief obtain the integration element
       *
       *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
       *  integration element \f$\mu(x)\f$ is given by
       *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
       *
       *  \param[in]  x  local coordinate to evaluate the integration element in
       *
       *  \returns the integration element \f$\mu(x)\f$.
       *
       *  \note For affine mappings, it is more efficient to call
       *        jacobianInverseTransposed before integrationElement, if both
       *        are required.
       */
      virtual FieldType integrationElement ( const LocalCoordinate &x ) const = 0;
      /** \brief obtain the volume of the mapping's image
       *
       *  \note The current implementation just returns
       *  \code
       *  integrationElement( baryCenter() ) * ReferenceElement::volume()
       *  \endcode
       *  which is wrong for n-linear surface maps and other nonlinear maps.
       */
      virtual FieldType volume () const = 0;

      /** \brief obtain the transposed of the Jacobian
       *
       *  \param[in]  x  local coordinate to evaluate Jacobian in
       *
       *  \returns a reference to the transposed of the Jacobian
       *
       *  \note The returned reference is reused on the next call to
       *        JacobianTransposed, destroying the previous value.
       */
      virtual const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &x ) const = 0;
      /** \brief obtain the transposed of the Jacobian's inverse
       *
       *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
       *  the Jacobian by \f$J(x)\f$, the following condition holds:
       *  \f[J^{-1}(x) J(x) = I.\f]
       */
      virtual const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &x ) const = 0;

    protected:
      using HybridMappingBase< dim, GeometryTraits >::trace;

    public:
      virtual This *clone () const = 0;
      virtual This *clone ( char *mappingStorage ) const = 0;

      template< int codim >
      typename Codim< codim >::Trace *trace ( unsigned int i, char *mappingStorage ) const
      {
        integral_constant< int, codim > codimVariable;
        return trace( codimVariable, i, mappingStorage );
      }

      const UserData &userData () const { return userData_; }
      UserData &userData () { return userData_; }

    private:
      UserData userData_;
    };



    // VirtualMappingBase
    // ------------------

    /** \cond */
    template< class Topology, class GeometryTraits, unsigned int codim = Topology::dimension >
    class VirtualMappingBase;

    template< class Topology, class GeometryTraits, unsigned int codim >
    class VirtualMappingBase
      : public VirtualMappingBase< Topology, GeometryTraits, codim-1 >,
        public virtual HybridMappingBase< Topology::dimension, GeometryTraits, codim >
    {
      typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits >
      VirtualMapping;

    protected:
      using VirtualMappingBase< Topology, GeometryTraits, codim-1 >::trace;

      virtual HybridMapping< Topology::dimension - codim, GeometryTraits > *
      trace ( integral_constant< int, codim >, unsigned int i, char *mappingStorage ) const
      {
        return static_cast< const VirtualMapping & >( *this ).template trace< codim >( i, mappingStorage );
      }
    };

    template< class Topology, class GeometryTraits >
    class VirtualMappingBase< Topology, GeometryTraits, 0 >
      : public virtual HybridMappingBase< Topology::dimension, GeometryTraits, 0 >
    {
      typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits >
      VirtualMapping;

    protected:
      virtual HybridMapping< Topology::dimension, GeometryTraits > *
      trace ( integral_constant< int, 0 >, unsigned int i, char *mappingStorage ) const
      {
        return static_cast< const VirtualMapping & >( *this ).template trace< 0 >( i, mappingStorage );
      }
    };
    /** \endcond */



    template< class Topology, class GeometryTraits >
    class VirtualMapping
      : public HybridMapping< Topology::dimension, GeometryTraits >,
        /** \cond */
        public VirtualMappingBase< Topology, GeometryTraits >
        /**\endcond*/
    {
      typedef HybridMapping< Topology::dimension, GeometryTraits > Base;
      typedef VirtualMapping< Topology, GeometryTraits > This;

      typedef typename Base::Traits Traits;

      typedef CachedMapping< Topology, GeometryTraits > Mapping;

    public:
      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      typedef typename Base::JacobianTransposed JacobianTransposed;
      typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;

      typedef typename Mapping::ReferenceElement ReferenceElement;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename TraceProvider< Topology, GeometryTraits, codim, true >::Trace Trace;
      };

      typedef typename GeometryTraits::Caching Caching;

      template< class CoordVector >
      explicit VirtualMapping ( const CoordVector &coordVector )
        : mapping_( coordVector )
      {}

      virtual bool affine () const { return mapping_.affine(); }
      virtual Dune::GeometryType type () const { return mapping_.type(); }

      virtual int numCorners () const { return mapping_.numCorners(); }
      virtual GlobalCoordinate corner ( int i ) const { return mapping_.corner( i ); }
      virtual GlobalCoordinate center () const { return mapping_.center(); }

      virtual GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_.global( local ); }
      virtual LocalCoordinate local ( const GlobalCoordinate &global ) const { return mapping_.local( global ); }

      virtual bool checkInside ( const LocalCoordinate &local ) const { return mapping_.checkInside( local ); }

      virtual FieldType integrationElement ( const LocalCoordinate &local ) const { return mapping_.integrationElement( local ); }
      virtual FieldType volume () const { return mapping_.volume(); }

      virtual const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_.jacobianTransposed( local ); }
      virtual const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_.jacobianInverseTransposed( local ); }

      virtual Base *clone () const { return new This( *this ); }
      virtual Base* clone ( char *mappingStorage ) const { return new( mappingStorage ) This( *this ); }

      template< int codim >
      typename Codim< codim >::Trace *trace ( unsigned int i, char *mappingStorage ) const
      {
        return TraceProvider< Topology, GeometryTraits, codim, true >::construct( mapping_.mapping(), i, mappingStorage );
      }

    protected:
      using VirtualMappingBase< Topology, GeometryTraits >::trace;

    private:
      Mapping mapping_;
    };



    // NonHybridMapping
    // ----------------

    /** \class   NonHybridMapping
     *  \ingroup GenericGeometry
     *  \brief   non-virtual geometric mapping
     *
     *  This is the user-visible class of the generic geometries if the
     *  topology type for each codimension is unique.
     */
    template< class Topology, class GeometryTraits >
    class NonHybridMapping
    {
      typedef NonHybridMapping< Topology, GeometryTraits > This;

    protected:
      typedef MappingTraits< typename GeometryTraits::CoordTraits, Topology::dimension, GeometryTraits::dimWorld > Traits;

      typedef CachedMapping< Topology, GeometryTraits > Mapping;

    public:
      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      typedef CachedJacobianTransposed< dimension, GeometryTraits > JacobianTransposed;
      typedef CachedJacobianInverseTransposed< dimension, GeometryTraits > JacobianInverseTransposed;

      typedef typename Mapping::ReferenceElement ReferenceElement;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename TraceProvider< Topology, GeometryTraits, codim, false >::Trace Trace;
      };

      typedef typename GeometryTraits::Caching Caching;
      typedef typename GeometryTraits::UserData UserData;

      template< class CoordVector >
      explicit NonHybridMapping ( const CoordVector &coordVector )
        : mapping_( coordVector )
      {}

      /** \brief is this mapping affine? */
      bool affine () const { return mapping_.affine(); }
      /** \brief obtain the name of the reference element */
      Dune::GeometryType type () const { return mapping_.type(); }

      /** \brief obtain number of corners of the corresponding reference element */
      int numCorners () const { return mapping_.numCorners(); }
      /** \brief obtain coordinates of the i-th corner */
      GlobalCoordinate corner ( int i ) const { return mapping_.corner( i ); }
      /** \brief obtain the centroid of the mapping's image */
      GlobalCoordinate center () const { return mapping_.center(); }

      /** \brief evaluate the mapping
       *
       *  \param[in]  local  local coordinate to map
       *
       *  \returns corresponding global coordinate
       */
      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_.global( local ); }
      /** \brief evaluate the inverse mapping
       *
       *  \param[in]  global  global coorindate to map
       *
       *  \return corresponding local coordinate
       *
       *  \note The returned local coordinate y minimizes
       *  \code
       *  (global( x ) - y).two_norm()
       *  \endcode
       */
      LocalCoordinate local ( const GlobalCoordinate &global ) const { return mapping_.local( global ); }

      /** \brief check whether a point lies within the reference element
       *
       *  \param[in]  local  local coordinate of point to check
       *
       *  \note Historically, this method was part of the geometry interface.
       *        It is still required for the GenericReferenceElement.
       */
      bool checkInside ( const LocalCoordinate &local ) const { return mapping_.checkInside( local ); }

      /** \brief obtain the integration element
       *
       *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
       *  integration element \f$\mu(x)\f$ is given by
       *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
       *
       *  \param[in]  local  local coordinate to evaluate the integration element in
       *
       *  \returns the integration element \f$\mu(x)\f$.
       *
       *  \note For affine mappings, it is more efficient to call
       *        jacobianInverseTransposed before integrationElement, if both
       *        are required.
       */
      FieldType integrationElement ( const LocalCoordinate &local ) const { return mapping_.integrationElement( local ); }
      /** \brief obtain the volume of the mapping's image
       *
       *  \note The current implementation just returns
       *  \code
       *  integrationElement( baryCenter() ) * ReferenceElement::volume()
       *  \endcode
       *  which is wrong for n-linear surface maps and other nonlinear maps.
       */
      FieldType volume () const { return mapping_.volume(); }

      /** \brief obtain the transposed of the Jacobian
       *
       *  \param[in]  local  local coordinate to evaluate Jacobian in
       *
       *  \returns a reference to the transposed of the Jacobian
       *
       *  \note The returned reference is reused on the next call to
       *        JacobianTransposed, destroying the previous value.
       */
      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_.jacobianTransposed( local ); }
      /** \brief obtain the transposed of the Jacobian's inverse
       *
       *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
       *  the Jacobian by \f$J(x)\f$, the following condition holds:
       *  \f[J^{-1}(x) J(x) = I.\f]
       */
      const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_.jacobianInverseTransposed( local ); }

      This *clone () const { return new This( *this ); }
      This *clone ( char *mappingStorage ) const { return new( mappingStorage ) This( *this ); }

      template< int codim >
      typename Codim< codim >::Trace *trace ( unsigned int i, char *mappingStorage ) const
      {
        return TraceProvider< Topology, GeometryTraits, codim, false >::construct( mapping_.mapping(), i, mappingStorage );
      }

      const UserData &userData () const { return userData_; }
      UserData &userData () { return userData_; }

    private:
      UserData userData_;
      Mapping mapping_;
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_HYBRIDMAPPING_HH
