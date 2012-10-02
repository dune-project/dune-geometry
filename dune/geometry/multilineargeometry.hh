// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH
#define DUNE_GEOMETRY_MULTILINEARGEOMETRY_HH

#include <cassert>
#include <limits>
#include <vector>

//#include <dune/common/array.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

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



  // MultiLinearGeometryTraits
  // -------------------------

  template< class ct >
  struct MultiLinearGeometryTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

    static ct tolerance () { return 16 * std::numeric_limits< ct >::epsilon(); }

    template< int mydim, int cdim >
    struct CornerStorage
    {
      //typedef array< FieldVector< ct, cdim >, 1 << mydim > Type;
      typedef std::vector< FieldVector< ct, cdim > > Type;
    };

    /** \brief will there be only one geometry type for a dimension?
     *
     *  If there is only a single geometry type for a certain dimension,
     *  <em>hasSingleGeometryType::v</em> can be set to true.
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

    struct UserData {};
  };



  // MultiLinearGeometry
  // -------------------

  template< class ct, int mydim, int cdim, class Traits = MultiLinearGeometryTraits< ct > >
  class MultiLinearGeometry
  {
    typedef MultiLinearGeometry< ct, mydim, cdim, Traits > This;

  public:
    typedef ct ctype;

    static const int mydimension= mydim;
    static const int coorddimension = cdim;

    typedef typename Traits::UserData UserData;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    class JacobianInverseTransposed;

    // for compatibility, export the type JacobianInverseTransposed as Jacobian
    typedef JacobianInverseTransposed Jacobian;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

  private:
    static const bool hasSingleGeometryType = Traits::template hasSingleGeometryType< mydimension >::v;

  protected:
    typedef typename Traits::MatrixHelper MatrixHelper;
    typedef typename SelectType< hasSingleGeometryType, integral_constant< unsigned int, Traits::template hasSingleGeometryType< mydimension >::topologyId >, unsigned int >::Type TopologyId;

    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

  private:
    struct Storage
      : public UserData
    {
      template< class CornerStorage >
      Storage ( const ReferenceElement &refEl, const CornerStorage &cornerStorage,
                const UserData &userData )
        : UserData( userData ),
          refElement( &refEl ),
          corners( cornerStorage )
      {}

      const ReferenceElement *refElement;
      typename Traits::template CornerStorage< mydimension, coorddimension >::Type corners;
    };

    typedef typename Traits::template CornerStorage< mydimension, coorddimension >::Type::const_iterator CornerIterator;

  public:
    template< class CornerStorage >
    MultiLinearGeometry ( const ReferenceElement &refElement, const CornerStorage &cornerStorage,
                          const UserData &userData = UserData() )
      : storage_( refElement, cornerStorage, userData )
    {}

    template< class CornerStorage >
    MultiLinearGeometry ( Dune::GeometryType gt, const CornerStorage &cornerStorage,
                          const UserData &userData = UserData() )
      : storage_( ReferenceElements::general( gt ), cornerStorage, userData )
    {}

    /** \brief is this mapping affine? */
    bool affine () const
    {
      CornerIterator cit = storage().corners.begin();
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
      return storage().corners[ i ];
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
      CornerIterator cit = storage().corners.begin();
      GlobalCoordinate y;
      global< false >( topologyId(), integral_constant< int, mydimension >(), cit, ctype( 1 ), local, ctype( 1 ), y );
      return y;
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
      const ctype tolerance = Traits::tolerance();
      LocalCoordinate x = refElement().position( 0, 0 );
      LocalCoordinate dx;
      do
      {
        // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
        const GlobalCoordinate dglobal = (*this).global( x ) - global;
        MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( x ), dglobal, dx );
        x -= dx;
      } while( dx.two_norm2() > tolerance*tolerance );
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
     *  integrationElement( baryCenter() ) * ReferenceElement::volume()
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
      CornerIterator cit = storage().corners.begin();
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

    const UserData &userData () const { return storage_; }
    UserData &userData () { return storage_; }

  protected:
    const Storage &storage () const { return storage_; }
    Storage &storage () { return storage_; }

    const ReferenceElement &refElement () const { return *storage().refElement; }

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

    template< bool add, int dim >
    static void jacobianTransposed ( TopologyId topologyId, integral_constant< int, dim >,
                                     CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                                     const ctype &rf, JacobianTransposed &jt );
    template< bool add >
    static void jacobianTransposed ( TopologyId topologyId, integral_constant< int, 0 >,
                                     CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                                     const ctype &rf, JacobianTransposed &jt );

    template< int dim >
    static bool affine ( TopologyId topologyId, integral_constant< int, dim >, CornerIterator &cit, JacobianTransposed &jt );
    static bool affine ( TopologyId topologyId, integral_constant< int, 0 >, CornerIterator &cit, JacobianTransposed &jt );

  protected:
    TopologyId topologyId ( integral_constant< bool, true > ) const { return TopologyId(); }
    unsigned int topologyId ( integral_constant< bool, false > ) const { return refElement().type().id(); }

    mutable JacobianTransposed jacobianTransposed_;
    mutable JacobianInverseTransposed jacobianInverseTransposed_;

  private:
    Storage storage_;
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



  // CachedMultiLinearGeometry
  // -------------------------

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
    typedef typename Base::UserData UserData;

    typedef typename Base::ctype ctype;

    using Base::mydimension;
    using Base::coorddimension;

    typedef typename Base::LocalCoordinate LocalCoordinate;
    typedef typename Base::GlobalCoordinate GlobalCoordinate;

    typedef typename Base::JacobianTransposed JacobianTransposed;
    typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;

    template< class CornerStorage >
    CachedMultiLinearGeometry ( const ReferenceElement &refElement, const CornerStorage &cornerStorage,
                                const UserData &userData = UserData() )
      : Base( refElement, cornerStorage, userData ),
        affine_( Base::affine() ),
        jacobianInverseTransposedComputed_( false ),
        integrationElementComputed_( false )
    {}

    template< class CornerStorage >
    CachedMultiLinearGeometry ( Dune::GeometryType gt, const CornerStorage &cornerStorage,
                                const UserData &userData = UserData() )
      : Base( gt, cornerStorage, userData ),
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
  template< bool add, int dim >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( TopologyId topologyId, integral_constant< int, dim >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, JacobianTransposed &jt )
  {
    const ctype xn = df*x[ dim-1 ];
    const ctype cxn = ctype( 1 ) - xn;
    assert( (xn > -Traits::tolerance()) && (cxn > -Traits::tolerance()) );

    if( GenericGeometry::isPrism( topologyId, mydimension, mydimension-dim ) )
    {
      CornerIterator cit2 = cit;
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
      CornerIterator cit2 = cit;
      // apply Jacobian for bottom (with argument x/(1-xn))
      jacobianTransposed< add >( topologyId, integral_constant< int, dim-1 >(), cit2, df/cxn, x, rf, jt );
      // compute last row
      global< add >( topologyId, integral_constant< int, dim-1 >(), cit, df/cxn, x, -rf, jt[ dim-1 ] );
      jt[ dim-1 ].axpy( rf, *cit );
      ++cit;
      for( int j = 0; j < dim-1; ++j )
        jt[ dim-1 ].axpy( rf*df*x[ j ], jt[ j ] );
    }
  }

  template< class ct, int mydim, int cdim, class Traits >
  template< bool add >
  inline void MultiLinearGeometry< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( TopologyId topologyId, integral_constant< int, 0 >,
                         CornerIterator &cit, const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, JacobianTransposed &jt )
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
