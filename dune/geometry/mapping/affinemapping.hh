// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_AFFINEMAPPING_HH
#define DUNE_GEOMETRY_AFFINEMAPPING_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElement;

  template< class ctype, int dim >
  class ReferenceElements;



  // AffineMappingTraits
  // -------------------

  template< class ct >
  struct AffineMappingTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

    struct UserData {};
  };



  // AffineMapping
  // -------------

  template< class ct, int mydim, int cdim, class Traits = AffineMappingTraits< ct > >
  class AffineMapping
  {
    typedef AffineMapping< ct, mydim, cdim > This;

  public:
    typedef ct ctype;

    static const int mydimension= mydim;
    static const int coorddimension = cdim;

    typedef typename Traits::UserData UserData;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    // for compatibility, export the type JacobianInverseTransposed as Jacobian
    typedef JacobianInverseTransposed Jacobian;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

  private:
    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

    typedef typename Traits::MatrixHelper MatrixHelper;

    struct Storage
      : public UserData
    {
      Storage ( const ReferenceElement &refEl, const GlobalCoordinate &org,
                const JacobianTransposed &jt, const UserData &userData )
        : UserData( userData ),
          refElement( &refEl ),
          origin( org ),
          jacobianTransposed( jt )
      {
        integrationElement = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed, jacobianInverseTransposed );
      }

      const ReferenceElement *refElement;
      GlobalCoordinate origin;
      JacobianTransposed jacobianTransposed;
      JacobianInverseTransposed jacobianInverseTransposed;
      ctype integrationElement;
    };

    AffineMapping ( const ReferenceElement &refElement, const GlobalCoordinate &origin,
                    const JacobianTransposed &jt, const UserData &userData = UserData() )
      : storage_( refElement, origin, jt, userData )
    {}

    AffineMapping ( Dune::GeometryType gt, const GlobalCoordinate &origin,
                    const JacobianTransposed &jt, const UserData &userData = UserData() )
      : storage_( ReferenceElements::general( gt ), origin, jt, userData )
    {}

    /** \brief is this mapping affine? */
    bool affine () const { return true; }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const { return refElement().type(); }

    /** \brief obtain number of corners of the corresponding reference element */
    int numCorners () const { return refElement().size( mydimension ); }

    /** \brief obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      return global( refElement().position( i, mydimension ) );
    }

    /** \brief obtain the centroid of the mapping's image */
    GlobalCoordinate center () const { return global( refElement().position( 0, 0 ) ); }

    /** \brief evaluate the mapping
     *
     *  \param[in]  local  local coordinate to map
     *
     *  \returns corresponding global coordinate
     */
    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      GlobalCoordinate global( storage().origin );
      storage().jacobianTransposed.umtv( local, global );
      return global;
    }

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
    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      LocalCoordinate local;
      storage().jacobianInverseTransposed.mtv( global - storage().origin, local );
      return local;
    }

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
    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return storage().integrationElement;
    }

    /** \brief obtain the volume of the mapping's image
     *
     *  \note The current implementation just returns
     *  \code
     *  integrationElement( baryCenter() ) * ReferenceElement::volume()
     *  \endcode
     *  which is wrong for n-linear surface maps and other nonlinear maps.
     */
    ctype volume () const
    {
      return storage().integrationElement * refElement().volume();
    }

    /** \brief obtain the transposed of the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a reference to the transposed of the Jacobian
     *
     *  \note The returned reference is reused on the next call to
     *        JacobianTransposed, destroying the previous value.
     */
    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
    {
      return storage().jacobianTranposed;
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return storage().jacobianInverseTransposed;
    }

    const UserData &userData () const { return storage_; }
    UserData &userData () { return storage_; }

  protected:
    const Storage &storage () const { return storage_; }
    Storage &storage () { return storage_; }

    const ReferenceElement &refElement () const { return *storage().refElement; }

  private:
    Storage storage_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_AFFINEMAPPING_HH
