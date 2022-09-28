// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH
#define DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH

#include <cassert>
#include <functional>
#include <iterator>
#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

namespace Dune
{

  // MultiLinearGeometryTraits
  // -------------------------

  /** \brief default traits class for MultiLinearGeometry
   *
   *  The MultiLinearGeometry (and CachedMultiLinearGeometry) allow tweaking
   *  some implementation details through a traits class.
   *
   *  This structure provides the default values.
   *
   *  \tparam  ct  coordinate type
   */
  template< class ct >
  struct MultiLinearGeometryTraits
  {
    /** \brief helper structure containing some matrix routines
     *
     *  This helper allows exchanging the matrix inversion algorithms.
     *  It must provide the following static methods:
     *  \code
     *  template< int m, int n >
     *  static ctype sqrtDetAAT ( const FieldMatrix< ctype, m, n > &A );
     *
     *  template< int m, int n >
     *  static ctype rightInvA ( const FieldMatrix< ctype, m, n > &A,
     *                           FieldMatrix< ctype, n, m > &ret );
     *
     *  template< int m, int n >
     *  static void xTRightInvA ( const FieldMatrix< ctype, m, n > &A,
     *                            const FieldVector< ctype, n > &x,
     *                            FieldVector< ctype, m > &y );
     *  \endcode
     */
    typedef Impl::FieldMatrixHelper< ct > MatrixHelper;

    /** \brief tolerance to numerical algorithms */
    static ct tolerance () { return ct( 16 ) * std::numeric_limits< ct >::epsilon(); }

    /** \brief template specifying the storage for the corners
     *
     *  Internally, the MultiLinearGeometry needs to store the corners of the
     *  geometry.
     *
     *  The corner storage may be chosen depending on geometry dimension and
     *  coordinate dimension. It is required to contain a type named Type, e.g.,
     *  \code
     *  template< int mydim, int cdim >
     *  struct CornerStorage
     *  {
     *    typedef std::vector< FieldVector< ctype, cdim > > Type;
     *  };
     *  \endcode
     *  By default, a std::vector of FieldVector is used.
     *
     *  Apart from being copy constructable and assignable, an \c const corner
     *  storage object \c corners must support the expressions \c
     *  begin(corners), \c end(corners), and subscription \c corners[i].  \c
     *  begin() and \c end() are looked up via ADL and in namespace \c std:
     *  \code
     *  using std::begin;
     *  using std::end;
     *  // it is a const_iterator over the corners in Dune-ordering
     *  auto it = begin(corners);
     *  FieldVector<ctype, cdim> c0 = *it;
     *  auto itend = end(corners);
     *  while(it != itend) {
     *    //...
     *  }
     *
     *  // elements must be accessible by subscription, indexed in
     *  // Dune-ordering
     *  FieldVector<ctype, cdim> c1 = corners[1];
     *  \endcode
     *  This means that all of the following qualify: \c
     *  FieldVector<ctype,cdim>[1<<mydim], \c
     *  std::array<FieldVector<ctype,cdim>,(1<<mydim)>, \c
     *  std::vector<FieldVector<ctype,cdim>>.
     *
     *  \note The expression \c end(corners) isn't actually used by the
     *        implementation currently, but we require it anyway so we can add
     *        runtime checks for the container size when we feel like it.
     *
     *  It is also possible to use a \c std::reference_wrapper of a suitable
     *  container as the type for the corner storage.  The implementation
     *  automatically calls \c corners.get() on internally stored \c
     *  std::reference_wrapper objects before applying \c begin(), \c end(),
     *  or subscription in that case.
     *
     *  \note Using \c std::reference_wrapper of some container as the corner
     *        storage means that the geometry has no control over the lifetime
     *        of or the access to that container.  When the lifetime of the
     *        container ends, or the container itself or its elements are
     *        modified, any geometry object that still references that
     *        container becomes invalid.  The only valid operation on invalid
     *        geometry objects are destruction and assignment from another
     *        geometry.  If invalidation happens concurrently with some
     *        operation (other than destruction or assignment) on the
     *        geometry, that is a race.
     *
     *  \tparam  mydim  geometry dimension
     *  \tparam  cdim   coordinate dimension
     */
    template< int mydim, int cdim >
    struct CornerStorage
    {
      typedef std::vector< FieldVector< ct, cdim > > Type;
    };

    /** \brief will there be only one geometry type for a dimension?
     *
     *  If there is only a single geometry type for a certain dimension,
     *  <em>hasSingleGeometryType::v</em> can be set to true.
     *  Supporting only one geometry type might yield a gain in performance.
     *
     *  If <em>hasSingleGeometryType::v</em> is set to true, an additional
     *  parameter <em>topologyId</em> is required.
     *  Here's an example:
     *  \code
     *  static const unsigned int topologyId = GeometryTypes::simplex(dim).id();
     *  \endcode
     */
    template< int dim >
    struct hasSingleGeometryType
    {
      static const bool v = false;
      static const unsigned int topologyId = ~0u;
    };
  };



  // MultiLinearGeometry
  // -------------------

  /** \brief generic geometry implementation based on corner coordinates
   *
   *  Based on the recursive definition of the reference elements, the
   *  MultiLinearGeometry provides a generic implementation of a geometry given
   *  the corner coordinates.
   *
   *  The geometric mapping is multilinear in the classical sense only in the
   *  case of cubes; for simplices it is linear.
   *  The name is still justified, because the mapping satisfies the important
   *  property of begin linear along edges.
   *
   *  \tparam  ct      coordinate type
   *  \tparam  mydim   geometry dimension
   *  \tparam  cdim    coordinate dimension
   *  \tparam  Traits  traits allowing to tweak some implementation details
   *                   (optional)
   *
   *  The requirements on the traits are documented along with their default,
   *  MultiLinearGeometryTraits.
   */
  template< class ct, int mydim, int cdim, class Traits = MultiLinearGeometryTraits< ct > >
  class MultiLinearGeometry
  {
    typedef MultiLinearGeometry< ct, mydim, cdim, Traits > This;

  public:
    //! coordinate type
    typedef ct ctype;

    //! geometry dimension
    static const int mydimension= mydim;
    //! coordinate dimension
    static const int coorddimension = cdim;

    //! type of local coordinates
    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    //! type of global coordinates
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
    //! type of volume
    typedef ctype Volume;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    //! type of jacobian inverse transposed
    class JacobianInverseTransposed;

    /** \brief Type for the Jacobian matrix */
    typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;

    /** \brief Type for the inverse Jacobian matrix */
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianInverse;

  protected:

    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

  public:

    //! type of reference element
    typedef typename ReferenceElements::ReferenceElement ReferenceElement;

  private:
    static const bool hasSingleGeometryType = Traits::template hasSingleGeometryType< mydimension >::v;

  protected:
    typedef typename Traits::MatrixHelper MatrixHelper;
    typedef typename std::conditional< hasSingleGeometryType, std::integral_constant< unsigned int, Traits::template hasSingleGeometryType< mydimension >::topologyId >, unsigned int >::type TopologyId;

  public:
    /** \brief constructor
     *
     *  \param[in]  refElement  reference element for the geometry
     *  \param[in]  corners     corners to store internally
     *
     *  \note The type of corners is actually a template argument.
     *        It is only required that the internal corner storage can be
     *        constructed from this object.
     */
    template< class Corners >
    MultiLinearGeometry ( const ReferenceElement &refElement,
                          const Corners &corners )
      : refElement_( refElement ),
        corners_( corners )
    {}

    /** \brief constructor
     *
     *  \param[in]  gt          geometry type
     *  \param[in]  corners     corners to store internally
     *
     *  \note The type of corners is actually a template argument.
     *        It is only required that the internal corner storage can be
     *        constructed from this object.
     */
    template< class Corners >
    MultiLinearGeometry ( Dune::GeometryType gt,
                          const Corners &corners )
      : refElement_( ReferenceElements::general( gt ) ),
        corners_( corners )
    {}

    /** \brief is this mapping affine? */
    bool affine () const
    {
      JacobianTransposed jt;
      return affine( jt );
    }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const { return GeometryType( toUnsignedInt(topologyId()), mydimension ); }

    /** \brief obtain number of corners of the corresponding reference element */
    int corners () const { return refElement().size( mydimension ); }

    /** \brief obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      assert( (i >= 0) && (i < corners()) );
      return std::cref(corners_).get()[ i ];
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
      using std::begin;

      auto cit = begin(std::cref(corners_).get());
      GlobalCoordinate y;
      global< false >( topologyId(), std::integral_constant< int, mydimension >(), cit, ctype( 1 ), local, ctype( 1 ), y );
      return y;
    }

    /** \brief evaluate the inverse mapping
     *
     *  \param[in] globalCoord global coordinate to map
     *
     *  \return corresponding local coordinate
     *
     *  \note For given global coordinate y the returned local coordinate x that minimizes
     *  the following function over the local coordinate space spanned by the reference element.
     *  \code
     *  (global( x ) - y).two_norm()
     *  \endcode
     */
    LocalCoordinate local ( const GlobalCoordinate &globalCoord ) const
    {
      const ctype tolerance = Traits::tolerance();
      LocalCoordinate x = refElement().position( 0, 0 );
      LocalCoordinate dx;
      const bool affineMapping = this->affine();
      do
      {
        // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
        const GlobalCoordinate dglobal = (*this).global( x ) - globalCoord;
        const bool invertible =
          MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( x ), dglobal, dx );
        if( ! invertible )
          return LocalCoordinate( std::numeric_limits< ctype > :: max() );

        // update x with correction
        x -= dx;

        // for affine mappings only one iteration is needed
        if ( affineMapping ) break;
      } while( dx.two_norm2() > tolerance );
      return x;
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
    Volume integrationElement ( const LocalCoordinate &local ) const
    {
      return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jacobianTransposed( local ) );
    }

    /** \brief obtain the volume of the mapping's image
     *
     *  \note The current implementation just returns
     *  \code
     *  integrationElement( refElement().position( 0, 0 ) ) * refElement().volume()
     *  \endcode
     *  which is wrong for n-linear surface maps and other nonlinear maps.
     */
    Volume volume () const
    {
      return integrationElement( refElement().position( 0, 0 ) ) * refElement().volume();
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
    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
      using std::begin;

      JacobianTransposed jt;
      auto cit = begin(std::cref(corners_).get());
      jacobianTransposed< false >( topologyId(), std::integral_constant< int, mydimension >(), cit, ctype( 1 ), local, ctype( 1 ), jt );
      return jt;
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const;

    friend ReferenceElement referenceElement ( const MultiLinearGeometry &geometry )
    {
      return geometry.refElement();
    }


    /** \brief Obtain the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a copy of the transposed of the Jacobian
     */
    Jacobian jacobian (const LocalCoordinate &local) const
    {
      return jacobianTransposed(local).transposed();
    }

    /** \brief Obtain the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    JacobianInverse jacobianInverse (const LocalCoordinate &local) const
    {
      return jacobianInverseTransposed(local).transposed();
    }

  protected:

    ReferenceElement refElement () const
    {
      return refElement_;
    }

    TopologyId topologyId () const
    {
      return topologyId( std::integral_constant< bool, hasSingleGeometryType >() );
    }

    template< bool add, int dim, class CornerIterator >
    static void global ( TopologyId topologyId, std::integral_constant< int, dim >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, GlobalCoordinate &y );
    template< bool add, class CornerIterator >
    static void global ( TopologyId topologyId, std::integral_constant< int, 0 >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, GlobalCoordinate &y );

    template< bool add, int rows, int dim, class CornerIterator >
    static void jacobianTransposed ( TopologyId topologyId, std::integral_constant< int, dim >,
                                     CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                                     const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt );
    template< bool add, int rows, class CornerIterator >
    static void jacobianTransposed ( TopologyId topologyId, std::integral_constant< int, 0 >,
                                     CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                                     const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt );

    template< int dim, class CornerIterator >
    static bool affine ( TopologyId topologyId, std::integral_constant< int, dim >, CornerIterator &cit, JacobianTransposed &jt );
    template< class CornerIterator >
    static bool affine ( TopologyId topologyId, std::integral_constant< int, 0 >, CornerIterator &cit, JacobianTransposed &jt );

    bool affine ( JacobianTransposed &jacobianT ) const
    {
      using std::begin;

      auto cit = begin(std::cref(corners_).get());
      return affine( topologyId(), std::integral_constant< int, mydimension >(), cit, jacobianT );
    }

  private:
    // The following methods are needed to convert the return type of topologyId to
    // unsigned int with g++-4.4. It has problems casting integral_constant to the
    // integral type.
    static unsigned int toUnsignedInt(unsigned int i) { return i; }
    template<unsigned int v>
    static unsigned int toUnsignedInt(std::integral_constant<unsigned int,v> ) { return v; }
    TopologyId topologyId ( std::integral_constant< bool, true > ) const { return TopologyId(); }
    unsigned int topologyId ( std::integral_constant< bool, false > ) const { return refElement().type().id(); }

    ReferenceElement refElement_;
    typename Traits::template CornerStorage< mydimension, coorddimension >::Type corners_;
  };



  // MultiLinearGeometry::JacobianInverseTransposed
  // ----------------------------------------------

  template< class ct, int mydim, int cdim, class Traits >
  class MultiLinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed
    : public FieldMatrix< ctype, coorddimension, mydimension >
  {
    typedef FieldMatrix< ctype, coorddimension, mydimension > Base;

  public:
    void setup ( const JacobianTransposed &jt )
    {
      detInv_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jt, static_cast< Base & >( *this ) );
    }

    void setupDeterminant ( const JacobianTransposed &jt )
    {
      detInv_ = MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jt );
    }

    ctype det () const { return ctype( 1 ) / detInv_; }
    ctype detInv () const { return detInv_; }

  private:
    ctype detInv_;
  };



  /** \brief Implement a MultiLinearGeometry with additional caching
   *
   * This class implements the same interface and functionality as MultiLinearGeometry.
   * However, it additionally implements caching for various results.
   *
   *  \tparam  ct      coordinate type
   *  \tparam  mydim   geometry dimension
   *  \tparam  cdim    coordinate dimension
   *  \tparam  Traits  traits allowing to tweak some implementation details
   *                   (optional)
   *
   */
  template< class ct, int mydim, int cdim, class Traits = MultiLinearGeometryTraits< ct > >
  class CachedMultiLinearGeometry
    : public MultiLinearGeometry< ct, mydim, cdim, Traits >
  {
    typedef CachedMultiLinearGeometry< ct, mydim, cdim, Traits > This;
    typedef MultiLinearGeometry< ct, mydim, cdim, Traits > Base;

  protected:
    typedef typename Base::MatrixHelper MatrixHelper;

  public:
    typedef typename Base::ReferenceElement ReferenceElement;

    typedef typename Base::ctype ctype;

    using Base::mydimension;
    using Base::coorddimension;

    typedef typename Base::LocalCoordinate LocalCoordinate;
    typedef typename Base::GlobalCoordinate GlobalCoordinate;
    typedef typename Base::Volume Volume;

    typedef typename Base::JacobianTransposed JacobianTransposed;
    typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;
    typedef typename Base::Jacobian Jacobian;
    typedef typename Base::JacobianInverse JacobianInverse;

    template< class CornerStorage >
    CachedMultiLinearGeometry ( const ReferenceElement &referenceElement, const CornerStorage &cornerStorage )
      : Base( referenceElement, cornerStorage ),
        affine_( Base::affine( jacobianTransposed_ ) ),
        jacobianInverseTransposedComputed_( false ),
        integrationElementComputed_( false )
    {}

    template< class CornerStorage >
    CachedMultiLinearGeometry ( Dune::GeometryType gt, const CornerStorage &cornerStorage )
      : Base( gt, cornerStorage ),
        affine_( Base::affine( jacobianTransposed_ ) ),
        jacobianInverseTransposedComputed_( false ),
        integrationElementComputed_( false )
    {}

    /** \brief is this mapping affine? */
    bool affine () const { return affine_; }

    using Base::corner;

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
      if( affine() )
      {
        GlobalCoordinate global( corner( 0 ) );
        jacobianTransposed_.umtv( local, global );
        return global;
      }
      else
        return Base::global( local );
    }

    /** \brief evaluate the inverse mapping
     *
     *  \param[in]  global  global coordinate to map
     *
     *  \return corresponding local coordinate
     *
     *  \note For given global coordinate y the returned local coordinate x that minimizes
     *  the following function over the local coordinate space spanned by the reference element.
     *  \code
     *  (global( x ) - y).two_norm()
     *  \endcode
     */
    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      if( affine() )
      {
        LocalCoordinate local;
        if( jacobianInverseTransposedComputed_ )
          jacobianInverseTransposed_.mtv( global - corner( 0 ), local );
        else
          MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed_, global - corner( 0 ), local );
        return local;
      }
      else
        return Base::local( global );
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
      if( affine() )
      {
        if( !integrationElementComputed_ )
        {
          jacobianInverseTransposed_.setupDeterminant( jacobianTransposed_ );
          integrationElementComputed_ = true;
        }
        return jacobianInverseTransposed_.detInv();
      }
      else
        return Base::integrationElement( local );
    }

    /** \brief obtain the volume of the mapping's image */
    Volume volume () const
    {
      if( affine() )
        return integrationElement( refElement().position( 0, 0 ) ) * refElement().volume();
      else
        return Base::volume();
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
    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
      if( affine() )
        return jacobianTransposed_;
      else
        return Base::jacobianTransposed( local );
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      if( affine() )
      {
        if( !jacobianInverseTransposedComputed_ )
        {
          jacobianInverseTransposed_.setup( jacobianTransposed_ );
          jacobianInverseTransposedComputed_ = true;
          integrationElementComputed_ = true;
        }
        return jacobianInverseTransposed_;
      }
      else
        return Base::jacobianInverseTransposed( local );
    }

    /** \brief Obtain the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a copy of the transposed of the Jacobian
     */
    Jacobian jacobian (const LocalCoordinate &local) const
    {
      return jacobianTransposed(local).transposed();
    }

    /** \brief Obtain the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    JacobianInverse jacobianInverse (const LocalCoordinate &local) const
    {
      return jacobianInverseTransposed(local).transposed();
    }

  protected:
    using Base::refElement;

  private:
    mutable JacobianTransposed jacobianTransposed_;
    mutable JacobianInverseTransposed jacobianInverseTransposed_;

    mutable bool affine_ : 1;

    mutable bool jacobianInverseTransposedComputed_ : 1;
    mutable bool integrationElementComputed_ : 1;
  };



  // Implementation of MultiLinearGeometry
  // -------------------------------------

  template< class ct, int mydim, int cdim, class Traits >
  inline typename MultiLinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed
  MultiLinearGeometry< ct, mydim, cdim, Traits >::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
    JacobianInverseTransposed jit;
    jit.setup( jacobianTransposed( local ) );
    return jit;
  }


  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, int dim, class CornerIterator >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::global ( TopologyId topologyId, std::integral_constant< int, dim >,
             CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
             const ctype &rf, GlobalCoordinate &y )
  {
    const ctype xn = df*x[ dim-1 ];
    const ctype cxn = ctype( 1 ) - xn;

    if( Impl::isPrism( toUnsignedInt(topologyId), mydimension, mydimension-dim ) )
    {
      // apply (1-xn) times mapping for bottom
      global< add >( topologyId, std::integral_constant< int, dim-1 >(), cit, df, x, rf*cxn, y );
      // apply xn times mapping for top
      global< true >( topologyId, std::integral_constant< int, dim-1 >(), cit, df, x, rf*xn, y );
    }
    else
    {
      assert( Impl::isPyramid( toUnsignedInt(topologyId), mydimension, mydimension-dim ) );
      // apply (1-xn) times mapping for bottom (with argument x/(1-xn))
      if( cxn > Traits::tolerance() || cxn < -Traits::tolerance() )
        global< add >( topologyId, std::integral_constant< int, dim-1 >(), cit, df/cxn, x, rf*cxn, y );
      else
        global< add >( topologyId, std::integral_constant< int, dim-1 >(), cit, df, x, ctype( 0 ), y );
      // apply xn times the tip
      y.axpy( rf*xn, *cit );
      ++cit;
    }
  }

  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, class CornerIterator >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::global ( TopologyId, std::integral_constant< int, 0 >,
             CornerIterator &cit, const ctype &, const LocalCoordinate &,
             const ctype &rf, GlobalCoordinate &y )
  {
    const GlobalCoordinate &origin = *cit;
    ++cit;
    for( int i = 0; i < coorddimension; ++i )
      y[ i ] = (add ? y[ i ] + rf*origin[ i ] : rf*origin[ i ]);
  }


  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, int rows, int dim, class CornerIterator >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( TopologyId topologyId, std::integral_constant< int, dim >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt )
  {
    assert( rows >= dim );

    const ctype xn = df*x[ dim-1 ];
    const ctype cxn = ctype( 1 ) - xn;

    auto cit2( cit );
    if( Impl::isPrism( toUnsignedInt(topologyId), mydimension, mydimension-dim ) )
    {
      // apply (1-xn) times Jacobian for bottom
      jacobianTransposed< add >( topologyId, std::integral_constant< int, dim-1 >(), cit2, df, x, rf*cxn, jt );
      // apply xn times Jacobian for top
      jacobianTransposed< true >( topologyId, std::integral_constant< int, dim-1 >(), cit2, df, x, rf*xn, jt );
      // compute last row as difference between top value and bottom value
      global< add >( topologyId, std::integral_constant< int, dim-1 >(), cit, df, x, -rf, jt[ dim-1 ] );
      global< true >( topologyId, std::integral_constant< int, dim-1 >(), cit, df, x, rf, jt[ dim-1 ] );
    }
    else
    {
      assert( Impl::isPyramid( toUnsignedInt(topologyId), mydimension, mydimension-dim ) );
      /*
       * In the pyramid case, we need a transformation Tb: B -> R^n for the
       * base B \subset R^{n-1}. The pyramid transformation is then defined as
       *   T: P \subset R^n  -> R^n
       *      (x, xn)        |-> (1-xn) Tb(x*) + xn t  (x \in R^{n-1}, xn \in R)
       * with the tip of the pyramid mapped to t and x* = x/(1-xn)
       * the projection of (x,xn) onto the base.
       *
       * For the Jacobi matrix DT we get
       *   DT = ( A | b )
       * with A = DTb(x*)   (n x n-1 matrix)
       *  and b = dT/dxn    (n-dim column vector).
       * Furthermore
       *   b = -Tb(x*) + t + \sum_i dTb/dx_i(x^*) x_i/(1-xn)
       *
       * Note that both A and b are not defined in the pyramid tip (x=0, xn=1)!
       * Indeed for B the unit square, Tb mapping B to the quadrilateral given
       * by the vertices (0,0,0), (2,0,0), (0,1,0), (1,1,0) and t=(0,0,1), we get
       *
       *   T(x,y,xn) = ( x(2-y/(1-xn)), y, xn )
       *               / 2-y/(1-xn)  -x   0 \
       *  DT(x,y,xn) = |    0         1   0 |
       *               \    0         0   1 /
       * which is not continuous for xn -> 1, choose for example
       *   x=0,    y=1-xn, xn -> 1   --> DT -> diag(1,1,1)
       *   x=1-xn, y=0,    xn -> 1   --> DT -> diag(2,1,1)
       *
       * However, for Tb affine-linear, Tb(y) = My + y0, DTb = M:
       *   A = M
       *   b = -M x* - y0 + t + \sum_i M_i x_i/(1-xn)
       *     = -M x* - y0 + t + M x*
       *     = -y0 + t
       * which is continuous for xn -> 1. Note that this b is also given by
       *   b = -Tb(0) + t + \sum_i dTb/dx_i(0) x_i/1
       * that is replacing x* by 1 and 1-xn by 1 in the formular above.
       *
       * For xn -> 1, we can thus set x*=0, "1-xn"=1 (or anything != 0) and get
       * the right result in case Tb is affine-linear.
       */

      /* The second case effectively results in x* = 0 */
      ctype dfcxn = (cxn > Traits::tolerance() || cxn < -Traits::tolerance()) ? ctype(df / cxn) : ctype(0);

      // initialize last row
      // b =  -Tb(x*)
      // (b = -Tb(0) = -y0 in case xn -> 1 and Tb affine-linear)
      global< add >( topologyId, std::integral_constant< int, dim-1 >(), cit, dfcxn, x, -rf, jt[ dim-1 ] );
      // b += t
      jt[ dim-1 ].axpy( rf, *cit );
      ++cit;
      // apply Jacobian for bottom (with argument x/(1-xn)) and correct last row
      if( add )
      {
        FieldMatrix< ctype, dim-1, coorddimension > jt2;
        // jt2 = dTb/dx_i(x*)
        jacobianTransposed< false >( topologyId, std::integral_constant< int, dim-1 >(), cit2, dfcxn, x, rf, jt2 );
        // A = dTb/dx_i(x*)                      (jt[j], j=0..dim-1)
        // b += \sum_i dTb/dx_i(x*) x_i/(1-xn)   (jt[dim-1])
        // (b += 0 in case xn -> 1)
        for( int j = 0; j < dim-1; ++j )
        {
          jt[ j ] += jt2[ j ];
          jt[ dim-1 ].axpy( dfcxn*x[ j ], jt2[ j ] );
        }
      }
      else
      {
        // jt = dTb/dx_i(x*)
        jacobianTransposed< false >( topologyId, std::integral_constant< int, dim-1 >(), cit2, dfcxn, x, rf, jt );
        // b += \sum_i dTb/dx_i(x*) x_i/(1-xn)
        for( int j = 0; j < dim-1; ++j )
          jt[ dim-1 ].axpy( dfcxn*x[ j ], jt[ j ] );
      }
    }
  }

  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, int rows, class CornerIterator >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( TopologyId, std::integral_constant< int, 0 >,
                         CornerIterator &cit, const ctype &, const LocalCoordinate &,
                         const ctype &, FieldMatrix< ctype, rows, cdim > & )
  {
    ++cit;
  }



  template< class ct, int mydim, int cdim, class Traits >
  template< int dim, class CornerIterator >
  inline bool MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::affine ( TopologyId topologyId, std::integral_constant< int, dim >, CornerIterator &cit, JacobianTransposed &jt )
  {
    const GlobalCoordinate &orgBottom = *cit;
    if( !affine( topologyId, std::integral_constant< int, dim-1 >(), cit, jt ) )
      return false;
    const GlobalCoordinate &orgTop = *cit;

    if( Impl::isPrism( toUnsignedInt(topologyId), mydimension, mydimension-dim ) )
    {
      JacobianTransposed jtTop;
      if( !affine( topologyId, std::integral_constant< int, dim-1 >(), cit, jtTop ) )
        return false;

      // check whether both jacobians are identical
      ctype norm( 0 );
      for( int i = 0; i < dim-1; ++i )
        norm += (jtTop[ i ] - jt[ i ]).two_norm2();
      if( norm >= Traits::tolerance() )
        return false;
    }
    else
      ++cit;
    jt[ dim-1 ] = orgTop - orgBottom;
    return true;
  }

  template< class ct, int mydim, int cdim, class Traits >
  template< class CornerIterator >
  inline bool MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::affine ( TopologyId, std::integral_constant< int, 0 >, CornerIterator &cit, JacobianTransposed & )
  {
    ++cit;
    return true;
  }

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH
