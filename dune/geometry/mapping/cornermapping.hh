// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_CORNERMAPPING_HH
#define DUNE_GEOMETRY_CORNERMAPPING_HH

#include <cassert>
#include <vector>

//#include <dune/common/array.hh>
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



  // CornerMappingTraits
  // -------------------

  template< class ct >
  struct CornerMappingTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

    static ct tolerance () { return 16 * std::numeric_limits< ct >::epsilon(); }

    template< int mydim, int cdim >
    struct CornerStorage
    {
      //typedef array< FieldVector< ct, cdim >, 1 << mydim > Type;
      typedef std::vector< FieldVector< ct, cdim > > Type;
    };

    struct UserData {};
  };



  // CornerMapping
  // -------------

  template< class ct, int mydim, int cdim, class Traits = CornerMappingTraits< ct > >
  class CornerMapping
  {
    typedef CornerMapping< ct, mydim, cdim, Traits > This;

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

  protected:
    typedef typename Traits::MatrixHelper MatrixHelper;

  private:
    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

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

  public:
    template< class CornerStorage >
    CornerMapping ( const ReferenceElement &refElement, const CornerStorage &cornerStorage,
                    const UserData &userData = UserData() )
      : storage_( refElement, cornerStorage, userData )
    {}

    template< class CornerStorage >
    CornerMapping ( Dune::GeometryType gt, const CornerStorage &cornerStorage,
                    const UserData &userData = UserData() )
      : storage_( ReferenceElements::general( gt ), cornerStorage, userData )
    {}

    /** \brief is this mapping affine? */
    bool affine () const
    {
      std::size_t offset = 0;
      return affine( type().id(), mydimension, offset, jacobianTransposed_ );
    }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const { return refElement().type(); }

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
      std::size_t offset = 0;
      GlobalCoordinate y;
      global< false >( type().id(), mydimension, offset, ctype( 1 ), local, ctype( 1 ), y );
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
      std::size_t offset = 0;
      jacobianTransposed< false >( type().id(), mydimension, offset, ctype( 1 ), local, ctype( 1 ), jacobianTransposed_ );
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

    template< bool add >
    void global ( unsigned int topologyId, int dim, std::size_t &offset,
                  const ctype &df, const LocalCoordinate &x,
                  const ctype &rf, GlobalCoordinate &y ) const;

    template< bool add >
    void jacobianTransposed ( unsigned int topologyId, int dim, std::size_t &offset,
                              const ctype &df, const LocalCoordinate &x,
                              const ctype &rf, JacobianTransposed &jt ) const;

    bool affine ( unsigned int topologyId, int dim, std::size_t &offset, JacobianTransposed &jt ) const;

  protected:
    mutable JacobianTransposed jacobianTransposed_;
    mutable JacobianInverseTransposed jacobianInverseTransposed_;

  private:
    Storage storage_;
  };



  // CornerMapping::JacobianInverseTransposed
  // ----------------------------------------

  template< class ct, int mydim, int cdim, class Traits >
  class CornerMapping< ct, mydim, cdim, Traits >::JacobianInverseTransposed
  {
    typedef FieldMatrix< ctype, coorddimension, mydimension > Matrix;

  public:
    static const int rows = Matrix::rows;
    static const int cols = Matrix::cols;

    void setup ( const JacobianTransposed &jt )
    {
      detInv_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jt, matrix_ );
    }

    void setupDeterminant ( const JacobianTransposed &jt )
    {
      detInv_ = MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jt );
    }

    operator const Matrix & () const { return matrix_; }

    template< class X, class Y >
    void mv ( const X &x, Y &y ) const { matrix_.mv( x, y ); }

    template< class X, class Y >
    void mtv ( const X &x, Y &y ) const { matrix_.mtv( x, y ); }

    template< class X, class Y >
    void umv ( const X &x, Y &y ) const { matrix_.umv( x, y ); }

    template< class X, class Y >
    void umtv ( const X &x, Y &y ) const { matrix_.umtv( x, y ); }

    template< class X, class Y >
    void mmv ( const X &x, Y &y ) const { matrix_.mmv( x, y ); }

    template< class X, class Y >
    void mmtv ( const X &x, Y &y ) const { matrix_.mmtv( x, y ); }

    ctype det () const { return ctype( 1 ) / detInv_; }
    ctype detInv () const { return detInv_; }

  private:
    Matrix matrix_;
    ctype detInv_;
  };



  // CachedCornerMapping
  // -------------------

  template< class ct, int mydim, int cdim, class Traits = CornerMappingTraits< ct > >
  class CachedCornerMapping
    : public CornerMapping< ct, mydim, cdim, Traits >
  {
    typedef CachedCornerMapping< ct, mydim, cdim, Traits > This;
    typedef CornerMapping< ct, mydim, cdim, Traits > Base;

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
    CachedCornerMapping ( const ReferenceElement &refElement, const CornerStorage &cornerStorage,
                          const UserData &userData = UserData() )
      : Base( refElement, cornerStorage, userData ),
        affine_( Base::affine() ),
        jacobianInverseTransposedComputed_( false ),
        integrationElementComputed_( false )
    {}

    template< class CornerStorage >
    CachedCornerMapping ( Dune::GeometryType gt, const CornerStorage &cornerStorage,
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



  // Implementation of CornerMapping
  // -------------------------------

  template< class ct, int mydim, int cdim, class Traits >
  inline const typename CornerMapping< ct, mydim, cdim, Traits >::JacobianInverseTransposed &
  CornerMapping< ct, mydim, cdim, Traits >
  ::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
    jacobianInverseTransposed_.setup( jacobianTransposed( local ) );
    return jacobianInverseTransposed_;
  }


  template< class ct, int mydim, int cdim, class Traits >
  template< bool add >
  inline void CornerMapping< ct, mydim, cdim, Traits >
  ::global ( unsigned int topologyId, int dim, std::size_t &offset,
             const ctype &df, const LocalCoordinate &x,
             const ctype &rf, GlobalCoordinate &y ) const
  {
    if( dim > 0 )
    {
      const unsigned int baseId = GenericGeometry::baseTopologyId( topologyId, dim );

      const ctype xn = df*x[ dim-1 ];
      const ctype cxn = ctype( 1 ) - xn;
      assert( (xn > -Traits::tolerance()) && (cxn > -Traits::tolerance()) );

      if( GenericGeometry::isPrism( topologyId, dim ) )
      {
        // apply (1-xn) times mapping for bottom
        global< add >( baseId, dim-1, offset, df, x, rf*cxn, y );
        // apply xn times mapping for top
        global< true >( baseId, dim-1, offset, df, x, rf*xn, y );
      }
      else
      {
        assert( GenericGeometry::isPyramid( topologyId, dim ) );
        // apply (1-xn) times mapping for bottom (with argument x/(1-xn))
        if( cxn > Traits::tolerance() )
          global< add >( baseId, dim-1, offset, df/cxn, x, rf*cxn, y );
        else
          global< add >( baseId, dim-1, offset, df, x, ctype( 0 ), y );
        // apply xn times the tip
        y.axpy( rf*xn, storage().corners[ offset++ ] );
      }
    }
    else
    {
      const GlobalCoordinate &origin = storage().corners[ offset++ ];
      for( int i = 0; i < coorddimension; ++i )
        y[ i ] = (add ? y[ i ] + rf*origin[ i ] : rf*origin[ i ]);
    }
  }


  template< class ct, int mydim, int cdim, class Traits >
  template< bool add >
  inline void CornerMapping< ct, mydim, cdim, Traits >
  ::jacobianTransposed ( unsigned int topologyId, int dim, std::size_t &offset,
                         const ctype &df, const LocalCoordinate &x,
                         const ctype &rf, JacobianTransposed &jt ) const
  {
    if( dim > 0 )
    {
      const unsigned int baseId = GenericGeometry::baseTopologyId( topologyId, dim );

      const ctype xn = df*x[ dim-1 ];
      const ctype cxn = ctype( 1 ) - xn;
      assert( (xn > -Traits::tolerance()) && (cxn > -Traits::tolerance()) );

      if( GenericGeometry::isPrism( topologyId, dim ) )
      {
        std::size_t it = offset;
        // apply (1-xn) times Jacobian for bottom
        jacobianTransposed< add >( baseId, dim-1, it, df, x, rf*cxn, jt );
        // apply xn times Jacobian for top
        jacobianTransposed< true >( baseId, dim-1, it, df, x, rf*xn, jt );
        // compute last row as difference between top value and bottom value
        global< add >( baseId, dim-1, offset, df, x, -rf, jt[ dim-1 ] );
        global< true >( baseId, dim-1, offset, df, x, rf, jt[ dim-1 ] );
      }
      else
      {
        assert( GenericGeometry::isPyramid( topologyId, dim ) );
        std::size_t it = offset;
        // apply Jacobian for bottom (with argument x/(1-xn))
        jacobianTransposed< add >( baseId, dim-1, it, df/cxn, x, rf, jt );
        // compute last row
        global< add >( baseId, dim-1, offset, df/cxn, x, -rf, jt[ dim-1 ] );
        jt[ dim-1 ].axpy( rf, storage().corners[ offset++ ] );
        for( int j = 0; j < dim-1; ++j )
          jt[ dim-1 ].axpy( rf*df*x[ j ], jt[ j ] );
      }
    }
    else
      ++offset;
  }


  template< class ct, int mydim, int cdim, class Traits >
  inline bool CornerMapping< ct, mydim, cdim, Traits >
  ::affine ( unsigned int topologyId, int dim, std::size_t &offset, JacobianTransposed &jt ) const
  {
    if( dim > 0 )
    {
      const unsigned int baseId = GenericGeometry::baseTopologyId( topologyId, dim );

      const GlobalCoordinate &orgBottom = storage().corners[ offset ];
      if( !affine( baseId, dim-1, offset, jt ) )
        return false;
      const GlobalCoordinate &orgTop = storage().corners[ offset ];

      if( GenericGeometry::isPrism( topologyId, dim ) )
      {
        JacobianTransposed jtTop;
        if( !affine( baseId, dim-1, offset, jtTop ) )
          return false;

        // check whether both jacobians are identical
        ctype norm( 0 );
        for( int i = 0; i < dim-1; ++i )
          norm += (jtTop[ i ] - jt[ i ]).two_norm2();
        if( norm >= Traits::tolerance() )
          return false;
      }
      jt[ dim-1 ] = orgTop - orgBottom;
    }
    else
      ++offset;
    return true;
  }

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_CORNERMAPPING_HH
