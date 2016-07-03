// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_AFFINEGEOMETRY_HH
#define DUNE_GEOMETRY_AFFINEGEOMETRY_HH

/** \file
 *  \brief An implementation of the Geometry interface for affine geometries
 *  \author Martin Nolte
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/implementation/matrixhelper.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElement;

  template< class ctype, int dim >
  struct ReferenceElements;



  /** \brief Implementation of the Geometry interface for affine geometries
   * \tparam ct Type used for coordinates
   * \tparam mydim Dimension of the geometry
   * \tparam cdim Dimension of the world space
   */
  template< class ct, int mydim, int cdim>
  class AffineGeometry
  {
  public:

    /** \brief Type used for coordinates */
    typedef ct ctype;

    /** \brief Dimension of the geometry */
    static const int mydimension= mydim;

    /** \brief Dimension of the world space */
    static const int coorddimension = cdim;

    /** \brief Type for local coordinate vector */
    typedef FieldVector< ctype, mydimension > LocalCoordinate;

    /** \brief Type for coordinate vector in world space */
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    /** \brief Type for the transposed Jacobian matrix */
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    /** \brief Type for the transposed inverse Jacobian matrix */
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

  private:
    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

    // Helper class to compute a matrix pseudo inverse
    typedef impl::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

  public:
    /** \brief Create affine geometry from reference element, one vertex, and the Jacobian matrix */
    AffineGeometry ( const ReferenceElement &refElement, const GlobalCoordinate &origin,
                     const JacobianTransposed &jt )
      : refElement_(&refElement), origin_(origin), jacobianTransposed_(jt)
    {
      integrationElement_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed_, jacobianInverseTransposed_ );
    }

    /** \brief Create affine geometry from GeometryType, one vertex, and the Jacobian matrix */
    AffineGeometry ( Dune::GeometryType gt, const GlobalCoordinate &origin,
                     const JacobianTransposed &jt )
      : refElement_( &ReferenceElements::general( gt ) ), origin_(origin), jacobianTransposed_( jt )
    {
      integrationElement_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed_, jacobianInverseTransposed_ );
    }

    /** \brief Create affine geometry from reference element and a vector of vertex coordinates */
    template< class CoordVector >
    AffineGeometry ( const ReferenceElement &refElement, const CoordVector &coordVector )
      : refElement_(&refElement), origin_(coordVector[0])
    {
      for( int i = 0; i < mydimension; ++i )
        jacobianTransposed_[ i ] = coordVector[ i+1 ] - origin_;
      integrationElement_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed_, jacobianInverseTransposed_ );
    }

    /** \brief Create affine geometry from GeometryType and a vector of vertex coordinates */
    template< class CoordVector >
    AffineGeometry ( Dune::GeometryType gt, const CoordVector &coordVector )
      : refElement_(&ReferenceElements::general( gt )), origin_(coordVector[0] )
    {
      for( int i = 0; i < mydimension; ++i )
        jacobianTransposed_[ i ] = coordVector[ i+1 ] - origin_;
      integrationElement_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed_, jacobianInverseTransposed_ );
    }

    /** \brief Always true: this is an affine geometry */
    bool affine () const { return true; }

    /** \brief Obtain the type of the reference element */
    Dune::GeometryType type () const { return refElement_->type(); }

    /** \brief Obtain number of corners of the corresponding reference element */
    int corners () const { return refElement_->size( mydimension ); }

    /** \brief Obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      return global( refElement_->position( i, mydimension ) );
    }

    /** \brief Obtain the centroid of the mapping's image */
    GlobalCoordinate center () const { return global( refElement_->position( 0, 0 ) ); }

    /** \brief Evaluate the mapping
     *
     *  \param[in]  local  local coordinate to map
     *
     *  \returns corresponding global coordinate
     */
    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      GlobalCoordinate global( origin_ );
      jacobianTransposed_.umtv( local, global );
      return global;
    }

    /** \brief Evaluate the inverse mapping
     *
     *  \param[in]  global  global coordinate to map
     *
     *  \return corresponding local coordinate
     *
     *  The returned local coordinate y minimizes
     *  \code
     *  (global( y ) - x).two_norm()
     *  \endcode
     *  on the entire affine hull of the reference element.  This degenerates
     *  to the inverse map if the argument y is in the range of the map.
     */
    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      LocalCoordinate local;
      jacobianInverseTransposed_.mtv( global - origin_, local );
      return local;
    }

    /** \brief Obtain the integration element
     *
     *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
     *  integration element \f$\mu(x)\f$ is given by
     *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
     *
     *  \param[in]  local  local coordinate to evaluate the integration element in
     *
     *  \returns the integration element \f$\mu(x)\f$.
     */
    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      DUNE_UNUSED_PARAMETER(local);
      return integrationElement_;
    }

    /** \brief Obtain the volume of the element */
    ctype volume () const
    {
      return integrationElement_ * refElement_->volume();
    }

    /** \brief Obtain the transposed of the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a reference to the transposed of the Jacobian
     */
    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
    {
      DUNE_UNUSED_PARAMETER(local);
      return jacobianTransposed_;
    }

    /** \brief Obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      DUNE_UNUSED_PARAMETER(local);
      return jacobianInverseTransposed_;
    }

    friend const ReferenceElement &referenceElement ( const AffineGeometry &geometry ) { return *geometry.refElement_; }

  private:
    const ReferenceElement* refElement_;
    GlobalCoordinate origin_;
    JacobianTransposed jacobianTransposed_;
    JacobianInverseTransposed jacobianInverseTransposed_;
    ctype integrationElement_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_AFFINEGEOMETRY_HH
