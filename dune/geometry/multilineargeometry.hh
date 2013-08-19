// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH
#define DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH

#include <cassert>
#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
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
  struct ReferenceElements;



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
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

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
     *  Apart from being copy constructable and assignable, the corner storage
     *  must provide a constant input iterator, i.e., it must define a type
     *  const_iterator and a pair of constant begin / end methods.
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
     *  static const unsigned int topologyId = SimplexTopology< dim >::type::id;
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
   *
   *  As an additional feature, this class allows to attach arbitrary user data
   *  to each object.  This is used in GeometryGrid to implement reference counting.
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

    /** \brief type of user data

        For example, GeometryGrid uses this to implement reference counting.
     */
    //! type of local coordinates
    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    //! type of global coordinates
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    //! type of jacobian inverse transposed
    class JacobianInverseTransposed;

    /** \brief For backward-compatibility, export the type JacobianInverseTransposed as Jacobian
     *  \deprecated This typedef will be removed after the release of dune 2.3
     */
    typedef JacobianInverseTransposed Jacobian;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

  private:
    static const bool hasSingleGeometryType = Traits::template hasSingleGeometryType< mydimension >::v;

  protected:
    typedef typename Traits::MatrixHelper MatrixHelper;
    typedef typename conditional< hasSingleGeometryType, integral_constant< unsigned int, Traits::template hasSingleGeometryType< mydimension >::topologyId >, unsigned int >::type TopologyId;

    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

  private:
    typedef typename Traits::template CornerStorage< mydimension, coorddimension >::Type::const_iterator CornerIterator;

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
      : refElement_( &refElement ),
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
      : refElement_( &ReferenceElements::general( gt ) ),
        corners_( corners )
    {}

    /** \brief is this mapping affine? */
    bool affine () const
    {
      CornerIterator cit = corners_.begin();
      return affine( topologyId(), integral_constant< int, mydimension >(), cit, jacobianTransposed_ );
    }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const { return GeometryType( topologyId(), mydimension ); }

    /** \brief obtain number of corners of the corresponding reference element */
    int corners () const { return refElement().size( mydimension ); }

    /** \brief obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      assert( (i >= 0) && (i < corners()) );
      return corners_[ i ];
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
      CornerIterator cit = corners_.begin();
      GlobalCoordinate y;
      global< false >( topologyId(), integral_constant< int, mydimension >(), cit, ctype( 1 ), local, ctype( 1 ), y );
      return y;
    }

    /** \brief evaluate the inverse mapping
     *
     *  \param[in]  global  global coordinate to map
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
      const ctype tolerance = Traits::tolerance();
      LocalCoordinate x = refElement().position( 0, 0 );
      LocalCoordinate dx;
      do
      {
        // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
        const GlobalCoordinate dglobal = (*this).global( x ) - global;
        MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( x ), dglobal, dx );
        x -= dx;
        assert( refElement().checkInside( x ) );
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
    ctype integrationElement ( const LocalCoordinate &local ) const
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
    ctype volume () const
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
    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
    {
      CornerIterator cit = corners_.begin();
      jacobianTransposed< false >( topologyId(), integral_constant< int, mydimension >(), cit, ctype( 1 ), local, ctype( 1 ), jacobianTransposed_ );
      return jacobianTransposed_;
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const;

  protected:
    const ReferenceElement &refElement () const { return *refElement_; }

    TopologyId topologyId () const
    {
      return topologyId( integral_constant< bool, hasSingleGeometryType >() );
    }

    template< bool add, int dim >
    static void global ( TopologyId topologyId, integral_constant< int, dim >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, GlobalCoordinate &y );
    template< bool add >
    static void global ( TopologyId topologyId, integral_constant< int, 0 >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, GlobalCoordinate &y );

    template< bool add, int rows, int dim >
    static void jacobianTransposed ( TopologyId topologyId, integral_constant< int, dim >,
                                     CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                                     const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt );
    template< bool add, int rows >
    static void jacobianTransposed ( TopologyId topologyId, integral_constant< int, 0 >,
                                     CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                                     const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt );

    template< int dim >
    static bool affine ( TopologyId topologyId, integral_constant< int, dim >, CornerIterator &cit, JacobianTransposed &jt );
    static bool affine ( TopologyId topologyId, integral_constant< int, 0 >, CornerIterator &cit, JacobianTransposed &jt );

  protected:
    TopologyId topologyId ( integral_constant< bool, true > ) const { return TopologyId(); }
    unsigned int topologyId ( integral_constant< bool, false > ) const { return refElement().type().id(); }

    mutable JacobianTransposed jacobianTransposed_;
    mutable JacobianInverseTransposed jacobianInverseTransposed_;

  private:
    const ReferenceElement *refElement_;
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

    typedef typename Base::JacobianTransposed JacobianTransposed;
    typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;

    template< class CornerStorage >
    CachedMultiLinearGeometry ( const ReferenceElement &refElement, const CornerStorage &cornerStorage )
      : Base( refElement, cornerStorage ),
        affine_( Base::affine() ),
        jacobianInverseTransposedComputed_( false ),
        integrationElementComputed_( false )
    {}

    template< class CornerStorage >
    CachedMultiLinearGeometry ( Dune::GeometryType gt, const CornerStorage &cornerStorage )
      : Base( gt, cornerStorage ),
        affine_( Base::affine() ),
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
     *  \note The returned local coordinate y minimizes
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
    ctype volume () const
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
    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
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
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const
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

  protected:
    using Base::refElement;

    using Base::jacobianTransposed_;
    using Base::jacobianInverseTransposed_;

  private:
    mutable bool affine_ : 1;

    mutable bool jacobianInverseTransposedComputed_ : 1;
    mutable bool integrationElementComputed_ : 1;
  };



  // Implementation of MultiLinearGeometry
  // -------------------------------------

  template< class ct, int mydim, int cdim, class Traits >
  inline const typename MultiLinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed &
  MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
    jacobianInverseTransposed_.setup( jacobianTransposed( local ) );
    return jacobianInverseTransposed_;
  }


  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, int dim >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::global ( TopologyId topologyId, integral_constant< int, dim >,
             CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
             const ctype &rf, GlobalCoordinate &y )
  {
    const ctype xn = df*x[ dim-1 ];
    const ctype cxn = ctype( 1 ) - xn;
    assert( (xn > -Traits::tolerance()) && (cxn > -Traits::tolerance()) );

    if( GenericGeometry::isPrism( topologyId, mydimension, mydimension-dim ) )
    {
      // apply (1-xn) times mapping for bottom
      global< add >( topologyId, integral_constant< int, dim-1 >(), cit, df, x, rf*cxn, y );
      // apply xn times mapping for top
      global< true >( topologyId, integral_constant< int, dim-1 >(), cit, df, x, rf*xn, y );
    }
    else
    {
      assert( GenericGeometry::isPyramid( topologyId, mydimension, mydimension-dim ) );
      // apply (1-xn) times mapping for bottom (with argument x/(1-xn))
      if( cxn > Traits::tolerance() )
        global< add >( topologyId, integral_constant< int, dim-1 >(), cit, df/cxn, x, rf*cxn, y );
      else
        global< add >( topologyId, integral_constant< int, dim-1 >(), cit, df, x, ctype( 0 ), y );
      // apply xn times the tip
      y.axpy( rf*xn, *cit );
      ++cit;
    }
  }

  template< class ct, int mydim, int cdim, class Traits >
  template< bool add >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::global ( TopologyId topologyId, integral_constant< int, 0 >,
             CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
             const ctype &rf, GlobalCoordinate &y )
  {
    const GlobalCoordinate &origin = *cit;
    ++cit;
    for( int i = 0; i < coorddimension; ++i )
      y[ i ] = (add ? y[ i ] + rf*origin[ i ] : rf*origin[ i ]);
  }


  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, int rows, int dim >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( TopologyId topologyId, integral_constant< int, dim >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt )
  {
    assert( rows >= dim );

    const ctype xn = df*x[ dim-1 ];
    const ctype cxn = ctype( 1 ) - xn;
    assert( (xn > -Traits::tolerance()) && (cxn > -Traits::tolerance()) );

    CornerIterator cit2( cit );
    if( GenericGeometry::isPrism( topologyId, mydimension, mydimension-dim ) )
    {
      // apply (1-xn) times Jacobian for bottom
      jacobianTransposed< add >( topologyId, integral_constant< int, dim-1 >(), cit2, df, x, rf*cxn, jt );
      // apply xn times Jacobian for top
      jacobianTransposed< true >( topologyId, integral_constant< int, dim-1 >(), cit2, df, x, rf*xn, jt );
      // compute last row as difference between top value and bottom value
      global< add >( topologyId, integral_constant< int, dim-1 >(), cit, df, x, -rf, jt[ dim-1 ] );
      global< true >( topologyId, integral_constant< int, dim-1 >(), cit, df, x, rf, jt[ dim-1 ] );
    }
    else
    {
      assert( GenericGeometry::isPyramid( topologyId, mydimension, mydimension-dim ) );
      // initialize last row
      global< add >( topologyId, integral_constant< int, dim-1 >(), cit, df/cxn, x, -rf, jt[ dim-1 ] );
      jt[ dim-1 ].axpy( rf, *cit );
      ++cit;
      // apply Jacobian for bottom (with argument x/(1-xn)) and correct last row
      if( add )
      {
        FieldMatrix< ctype, dim-1, coorddimension > jt2;
        jacobianTransposed< false >( topologyId, integral_constant< int, dim-1 >(), cit2, df/cxn, x, rf, jt2 );
        for( int j = 0; j < dim-1; ++j )
        {
          jt[ j ] += jt2[ j ];
          jt[ dim-1 ].axpy( (df/cxn)*x[ j ], jt2[ j ] );
        }
      }
      else
      {
        jacobianTransposed< false >( topologyId, integral_constant< int, dim-1 >(), cit2, df/cxn, x, rf, jt );
        for( int j = 0; j < dim-1; ++j )
          jt[ dim-1 ].axpy( (df/cxn)*x[ j ], jt[ j ] );
      }
    }
  }

  template< class ct, int mydim, int cdim, class Traits >
  template< bool add, int rows >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( TopologyId topologyId, integral_constant< int, 0 >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, FieldMatrix< ctype, rows, cdim > &jt )
  {
    ++cit;
  }



  template< class ct, int mydim, int cdim, class Traits >
  template< int dim >
  inline bool MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::affine ( TopologyId topologyId, integral_constant< int, dim >, CornerIterator &cit, JacobianTransposed &jt )
  {
    const GlobalCoordinate &orgBottom = *cit;
    if( !affine( topologyId, integral_constant< int, dim-1 >(), cit, jt ) )
      return false;
    const GlobalCoordinate &orgTop = *cit;

    if( GenericGeometry::isPrism( topologyId, mydimension, mydimension-dim ) )
    {
      JacobianTransposed jtTop;
      if( !affine( topologyId, integral_constant< int, dim-1 >(), cit, jtTop ) )
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
  inline bool MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::affine ( TopologyId topologyId, integral_constant< int, 0 >, CornerIterator &cit, JacobianTransposed &jt )
  {
    ++cit;
    return true;
  }

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH
