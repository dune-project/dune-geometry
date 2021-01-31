// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GEOMETRY_QUADRATURERULES_HH
#define DUNE_GEOMETRY_QUADRATURERULES_HH

#include <algorithm>
#include <iostream>
#include <limits>
#include <mutex>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/stdthread.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

/**
   \file
   Interface for quadrature points and rules
 */

namespace Dune {

  /** \brief Exception thrown if a desired QuadratureRule is not available,
     because the requested order is to high
     \ingroup Quadrature
   */
  class QuadratureOrderOutOfRange : public NotImplemented {};

  /** \brief Single evaluation point in a quadrature rule
      \ingroup Quadrature
      \tparam ct Number type used for both coordinates and the weights
      \tparam dim Dimension of the integration domain
   */
  template<typename ct, int dim>
  class QuadraturePoint {
  public:
    /** \brief Dimension of the integration domain */
    enum { dimension = dim };

    /** \brief Number type used for coordinates and quadrature weights */
    typedef ct Field;

    /** \brief Type used for the position of a quadrature point */
    typedef Dune::FieldVector<ct,dim> Vector;

    //! set up quadrature of given order in d dimensions
    QuadraturePoint (const Vector& x, ct w) : local(x)
    {
      weight_ = w;
    }

    //! return local coordinates of integration point i
    const Vector& position () const
    {
      return local;
    }

    //! return weight associated with integration point i
    const ct &weight () const
    {
      return weight_;
    }

  protected:
    FieldVector<ct, dim> local;
    ct weight_;
  };

  /** \brief Defines an \p enum for currently available quadrature rules.
      \ingroup Quadrature
   */
  namespace QuadratureType {
    enum Enum {
      /** \brief Gauss-Legendre rules (default)
      *
      *  -1D: Gauss-Jacobi rule with parameters \f$\alpha = \beta =0 \f$, i.e. for integrals with a constant weight function.
      *       The quadrature points do not include interval endpoints.
      *       Polynomials of order 2n - 1 can be integrated exactly.
      *  -higher dimension: For the 2D/3D case efficient rules for certain geometries may be used if available.
      *                     Higher dimensional quadrature rules are constructed via \p TensorProductQuadratureRule.
      *                     In this case the 1D rules eventually need higher order to compensate occuring weight functions(i.e. simplices).
      */
      GaussLegendre = 0,

      /** \brief Gauss-Jacobi rules with \f$\alpha =1\f$
      *
      *  -1D Gauss-Jacobi rule with parameters \f$\alpha =1,\ \beta =0 \f$
      *  -Is used to construct efficient simplex quadrature rules of higher order
      */
      GaussJacobi_1_0 = 1,

      /** \brief Gauss-Legendre rules with \f$\alpha =2\f$
      *
      *  -1D Gauss-Jacobi rule with parameters \f$\alpha =2,\ \beta =0 \f$
      *  -Is used to construct efficient simplex quadrature rules of higher order
      */
      GaussJacobi_2_0 = 2,

      /** \brief Gauss-Legendre rules with \f$\alpha =n\f$
      *
      *  -1D: Gauss-Jacobi rule with parameters \f$\alpha = n,\ \beta =0 \f$
      *  -higher dimension: For the 2D/3D case efficient rules for certain geometries may be used if available.
      *                     Higher dimensional quadrature rules are constructed via \p TensorProductQuadratureRule.
      *                     In this case the 1D rules respect eventually occuring weight functions(i.e. simplices).
      *  -The rules for high dimension or order are computed at run time and only floating point number types are supported.(LAPACK is needed for this case)
      *  -Most efficient quadrature type for simplices.
      *
      *   \note For details please use the book "Approximate Calculation of Multiple Integrals" by A.H. Stroud published in 1971.
      */
      GaussJacobi_n_0 = 3,

      /** \brief Gauss-Lobatto rules
       *
       * 1D: Gauss-Lobatto rules for a constant weight function.
       * These are optimal rules under the constraint that both interval endpoints are quadrature points.
       * Polynomials of order 2n - 3 can be integrated exactly.
       */
      GaussLobatto = 4,

      /** \brief Gauss-Radau rules including the left endpoint
       *
       * 1D: Gauss-Radau rules for a constant weight function.
       * These are optimal rules under the constraint that the left endpoint of the integration interval is a quadrature point.
       * Polynomials of order 2n - 2 can be integrated exactly.
       */
      GaussRadauLeft = 5,

      /** \brief Gauss-Radau rules including the right endpoint
       *
       * 1D: Gauss-Radau rules for a constant weight function.
       * These are optimal rules under the constraint that the right endpoint of the integration interval is a quadrature point.
       * Polynomials of order 2n - 2 can be integrated exactly.
       * The right Gauss-Radau rules are the just the mirrored left Gauss-Radau rules.
       */
      GaussRadauRight = 6,
      size
    };
  }

  /** \brief Abstract base class for quadrature rules
      \ingroup Quadrature
   */
  template<typename ct, int dim>
  class QuadratureRule : public std::vector<QuadraturePoint<ct,dim> >
  {
  public:
    /** \brief Default constructor
     *
     * Create an invalid empty quadrature rule.  This must be initialized
     * later by copying another quadraturerule before it can be used.
     */
    QuadratureRule() : delivered_order(-1) {}

  protected:
    /** \brief Constructor for a given geometry type.  Leaves the quadrature order invalid  */
    QuadratureRule(GeometryType t) : geometry_type(t), delivered_order(-1) {}

    /** \brief Constructor for a given geometry type and a given quadrature order */
    QuadratureRule(GeometryType t, int order) : geometry_type(t), delivered_order(order) {}
  public:
    /** \brief The space dimension */
    enum { d=dim };

    /** \brief The type used for coordinates */
    typedef ct CoordType;

    //! return order
    virtual int order () const { return delivered_order; }

    //! return type of element
    virtual GeometryType type () const { return geometry_type; }
    virtual ~QuadratureRule(){}

    //! this container is always a const container,
    //! therefore iterator is the same as const_iterator
    typedef typename std::vector<QuadraturePoint<ct,dim> >::const_iterator iterator;

  protected:
    GeometryType geometry_type;
    int delivered_order;
  };

  // Forward declaration of the factory class,
  // needed internally by the QuadratureRules container class.
  template<typename ctype, int dim> class QuadratureRuleFactory;

  /** \brief A container for all quadrature rules of dimension <tt>dim</tt>
      \ingroup Quadrature
   */
  template<typename ctype, int dim>
  class QuadratureRules {

    /** \brief Internal short-hand notation for the type of quadrature rules this container contains */
    typedef Dune::QuadratureRule<ctype, dim> QuadratureRule;
    //! \brief a quadrature rule (for each quadrature order, geometry type,
    //!        and quadrature type)
    static void initQuadratureRule(QuadratureRule *qr, QuadratureType::Enum qt,
                                   const GeometryType &t, int p)
    {
      *qr = QuadratureRuleFactory<ctype,dim>::rule(t,p,qt);
    }

    typedef std::vector<std::pair<std::once_flag, QuadratureRule> >
      QuadratureOrderVector; // indexed by quadrature order
    //! \brief initialize the vector indexed by the quadrature order (for each
    //!        geometry type and quadrature type)
    static void initQuadratureOrderVector(QuadratureOrderVector *qov,
                                          QuadratureType::Enum qt,
                                          const GeometryType &t)
    {
      if(dim == 0)
        // we only need one quadrature rule for points, not maxint
        *qov = QuadratureOrderVector(1);
      else
        *qov = QuadratureOrderVector(QuadratureRuleFactory<ctype,dim>::maxOrder(t,qt)+1);
    }

    typedef std::vector<std::pair<std::once_flag, QuadratureOrderVector> >
      GeometryTypeVector; // indexed by geometry type
    //! \brief initialize the vector indexed by the geometry type (for each
    //!        quadrature type)
    static void initGeometryTypeVector(GeometryTypeVector *gtv)
    {
      *gtv = GeometryTypeVector(LocalGeometryTypeIndex::size(dim));
    }

    //! real rule creator
    DUNE_EXPORT const QuadratureRule& _rule(const GeometryType& t, int p, QuadratureType::Enum qt=QuadratureType::GaussLegendre)
    {
      assert(t.dim()==dim);

      DUNE_ASSERT_CALL_ONCE();

      static std::vector<std::pair< // indexed by quadrature type
        std::once_flag,
        GeometryTypeVector
        > > quadratureCache(QuadratureType::size);

      auto & quadratureTypeLevel = quadratureCache[qt];
      std::call_once(quadratureTypeLevel.first, initGeometryTypeVector,
                     &quadratureTypeLevel.second);

      auto & geometryTypeLevel =
        quadratureTypeLevel.second[LocalGeometryTypeIndex::index(t)];
      std::call_once(geometryTypeLevel.first, initQuadratureOrderVector,
                     &geometryTypeLevel.second, qt, t);

      // we only have one quadrature rule for points
      auto & quadratureOrderLevel = geometryTypeLevel.second[dim == 0 ? 0 : p];
      std::call_once(quadratureOrderLevel.first, initQuadratureRule,
                     &quadratureOrderLevel.second, qt, t, p);

      return quadratureOrderLevel.second;
    }
    //! singleton provider
    DUNE_EXPORT static QuadratureRules& instance()
    {
      static QuadratureRules instance;
      return instance;
    }
    //! private constructor
    QuadratureRules () {}
  public:
    //! maximum quadrature order for given geometry type and quadrature type
    static unsigned
    maxOrder(const GeometryType& t,
             QuadratureType::Enum qt=QuadratureType::GaussLegendre)
    {
      return QuadratureRuleFactory<ctype,dim>::maxOrder(t,qt);
    }

    //! select the appropriate QuadratureRule for GeometryType t and order p
    static const QuadratureRule& rule(const GeometryType& t, int p, QuadratureType::Enum qt=QuadratureType::GaussLegendre)
    {
      return instance()._rule(t,p,qt);
    }

    //! @copydoc rule
    static const QuadratureRule& rule(const GeometryType::BasicType t, int p, QuadratureType::Enum qt=QuadratureType::GaussLegendre)
    {
      GeometryType gt(t,dim);
      return instance()._rule(gt,p,qt);
    }
  };

} // end namespace Dune

#include "quadraturerules/pointquadrature.hh"

namespace Dune {

  //! \internal Helper template for the initialization of the quadrature rules
  template<typename ct, bool fundamental = std::numeric_limits<ct>::is_specialized>
  struct GaussQuadratureInitHelper;
  template<typename ct>
  struct GaussQuadratureInitHelper<ct, true> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };
  template<typename ct>
  struct GaussQuadratureInitHelper<ct, false> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };

  //! \brief Gauss quadrature rule in 1D
  template<typename ct>
  class GaussQuadratureRule1D :
    public QuadratureRule<ct,1>
  {
  public:
    // compile time parameters
    enum { dim=1 };
    enum { highest_order=61 };

    ~GaussQuadratureRule1D(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    GaussQuadratureRule1D (int p)
      : QuadratureRule<ct,1>(GeometryTypes::line)
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      GaussQuadratureInitHelper<ct>::init
        (p, _points, _weight, this->delivered_order);

      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };

  extern template GaussQuadratureRule1D<float>::GaussQuadratureRule1D(int);
  extern template GaussQuadratureRule1D<double>::GaussQuadratureRule1D(int);

} // namespace Dune

#define DUNE_INCLUDING_IMPLEMENTATION
#include "quadraturerules/gauss_imp.hh"

namespace Dune {

  //! \internal Helper template for the initialization of the quadrature rules
  template<typename ct,
      bool fundamental = std::numeric_limits<ct>::is_specialized>
  struct Jacobi1QuadratureInitHelper;
  template<typename ct>
  struct Jacobi1QuadratureInitHelper<ct, true> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };
  template<typename ct>
  struct Jacobi1QuadratureInitHelper<ct, false> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };

  /** \brief Jacobi-Gauss quadrature for alpha=1, beta=0
      \ingroup Quadrature
   */
  template<typename ct>
  class Jacobi1QuadratureRule1D :
    public QuadratureRule<ct,1>
  {
  public:
    /** \brief The space dimension */
    enum { dim=1 };

    /** \brief The highest quadrature order available */
    enum { highest_order=61 };

    ~Jacobi1QuadratureRule1D(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    Jacobi1QuadratureRule1D (int p)
      : QuadratureRule<ct,1>(GeometryTypes::line)
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      int deliveredOrder_;

      Jacobi1QuadratureInitHelper<ct>::init
        (p, _points, _weight, deliveredOrder_);
      this->delivered_order = deliveredOrder_;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };

#ifndef DOXYGEN
  extern template Jacobi1QuadratureRule1D<float>::Jacobi1QuadratureRule1D(int);
  extern template Jacobi1QuadratureRule1D<double>::Jacobi1QuadratureRule1D(int);
#endif // !DOXYGEN

} // namespace Dune

#define DUNE_INCLUDING_IMPLEMENTATION
#include "quadraturerules/jacobi_1_0_imp.hh"

namespace Dune {

  //! \internal Helper template for the initialization of the quadrature rules
  template<typename ct,
      bool fundamental = std::numeric_limits<ct>::is_specialized>
  struct Jacobi2QuadratureInitHelper;
  template<typename ct>
  struct Jacobi2QuadratureInitHelper<ct, true> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };
  template<typename ct>
  struct Jacobi2QuadratureInitHelper<ct, false> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };

  /** \brief Jacobi-Gauss quadrature for alpha=2, beta=0
      \ingroup Quadrature
   */
  template<typename ct>
  class Jacobi2QuadratureRule1D :
    public QuadratureRule<ct,1>
  {
  public:
    /** \brief The space dimension */
    enum { dim=1 };

    /** \brief The highest quadrature order available */
    enum { highest_order=61 };

    ~Jacobi2QuadratureRule1D(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    Jacobi2QuadratureRule1D (int p)
      : QuadratureRule<ct,1>(GeometryTypes::line)
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      int deliveredOrder_;

      Jacobi2QuadratureInitHelper<ct>::init
        (p, _points, _weight, deliveredOrder_);

      this->delivered_order = deliveredOrder_;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };

#ifndef DOXYGEN
  extern template Jacobi2QuadratureRule1D<float>::Jacobi2QuadratureRule1D(int);
  extern template Jacobi2QuadratureRule1D<double>::Jacobi2QuadratureRule1D(int);
#endif // !DOXYGEN

} // namespace Dune

#define DUNE_INCLUDING_IMPLEMENTATION
#include "quadraturerules/jacobi_2_0_imp.hh"

namespace Dune {

  //! \internal Helper template for the initialization of the quadrature rules
  template<typename ct,
      bool fundamental = std::numeric_limits<ct>::is_specialized>
  struct GaussLobattoQuadratureInitHelper;
  template<typename ct>
  struct GaussLobattoQuadratureInitHelper<ct, true> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };
  template<typename ct>
  struct GaussLobattoQuadratureInitHelper<ct, false> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };

  /** \brief Jacobi-Gauss quadrature for alpha=2, beta=0
      \ingroup Quadrature
   */
  template<typename ct>
  class GaussLobattoQuadratureRule1D :
    public QuadratureRule<ct,1>
  {
  public:
    /** \brief The space dimension */
    enum { dim=1 };

    /** \brief The highest quadrature order available */
    enum { highest_order=31 };

    ~GaussLobattoQuadratureRule1D(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    GaussLobattoQuadratureRule1D (int p)
      : QuadratureRule<ct,1>(GeometryTypes::line)
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      int deliveredOrder_;

      GaussLobattoQuadratureInitHelper<ct>::init
        (p, _points, _weight, deliveredOrder_);

      this->delivered_order = deliveredOrder_;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };

#ifndef DOXYGEN
  extern template GaussLobattoQuadratureRule1D<float>::GaussLobattoQuadratureRule1D(int);
  extern template GaussLobattoQuadratureRule1D<double>::GaussLobattoQuadratureRule1D(int);
#endif // !DOXYGEN

} // namespace Dune

#define DUNE_INCLUDING_IMPLEMENTATION
#include "quadraturerules/gausslobatto_imp.hh"





namespace Dune {

  //! \internal Helper template for the initialization of the quadrature rules
  template<typename ct,
      bool fundamental = std::numeric_limits<ct>::is_specialized>
  struct GaussRadauLeftQuadratureInitHelper;
  template<typename ct>
  struct GaussRadauLeftQuadratureInitHelper<ct, true> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };
  template<typename ct>
  struct GaussRadauLeftQuadratureInitHelper<ct, false> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };

  /** \brief Gauss-Radau quadrature, including the left endpoint of the interval
      \ingroup Quadrature
   */
  template<typename ct>
  class GaussRadauLeftQuadratureRule1D :
    public QuadratureRule<ct,1>
  {
  public:
    /** \brief The space dimension */
    enum { dim=1 };

    /** \brief The highest quadrature order available */
    enum { highest_order=30};

    ~GaussRadauLeftQuadratureRule1D(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    GaussRadauLeftQuadratureRule1D (int p)
      : QuadratureRule<ct,1>(GeometryTypes::line)
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      int deliveredOrder_;

      GaussRadauLeftQuadratureInitHelper<ct>::init
        (p, _points, _weight, deliveredOrder_);

      this->delivered_order = deliveredOrder_;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };

#ifndef DOXYGEN
  extern template GaussRadauLeftQuadratureRule1D<float>::GaussRadauLeftQuadratureRule1D(int);
  extern template GaussRadauLeftQuadratureRule1D<double>::GaussRadauLeftQuadratureRule1D(int);
#endif // !DOXYGEN

} // namespace Dune

#define DUNE_INCLUDING_IMPLEMENTATION
#include "quadraturerules/gaussradauleft_imp.hh"



namespace Dune {

  //! \internal Helper template for the initialization of the quadrature rules
  template<typename ct,
      bool fundamental = std::numeric_limits<ct>::is_specialized>
  struct GaussRadauRightQuadratureInitHelper;
  template<typename ct>
  struct GaussRadauRightQuadratureInitHelper<ct, true> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };
  template<typename ct>
  struct GaussRadauRightQuadratureInitHelper<ct, false> {
    static void init(int p,
                     std::vector< FieldVector<ct, 1> > & _points,
                     std::vector< ct > & _weight,
                     int & delivered_order);
  };

  /** \brief Gauss-Radau quadrature, including the left endpoint of the interval
      \ingroup Quadrature
   */
  template<typename ct>
  class GaussRadauRightQuadratureRule1D :
    public QuadratureRule<ct,1>
  {
  public:
    /** \brief The space dimension */
    enum { dim=1 };

    /** \brief The highest quadrature order available */
    enum { highest_order=30};

    ~GaussRadauRightQuadratureRule1D(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    GaussRadauRightQuadratureRule1D (int p)
      : QuadratureRule<ct,1>(GeometryTypes::line)
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      int deliveredOrder_;

      GaussRadauLeftQuadratureInitHelper<ct>::init
        (p, _points, _weight, deliveredOrder_);

      this->delivered_order = deliveredOrder_;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };

#ifndef DOXYGEN
  extern template GaussRadauRightQuadratureRule1D<float>::GaussRadauRightQuadratureRule1D(int);
  extern template GaussRadauRightQuadratureRule1D<double>::GaussRadauRightQuadratureRule1D(int);
#endif // !DOXYGEN

} // namespace Dune

#define DUNE_INCLUDING_IMPLEMENTATION
#include "quadraturerules/gaussradauright_imp.hh"




#include "quadraturerules/tensorproductquadrature.hh"

#include "quadraturerules/simplexquadrature.hh"

namespace Dune {

  /***********************************
   * quadrature for Prism
   **********************************/

  /** \todo Please doc me! */
  template<int dim>
  class PrismQuadraturePoints;

  /** \todo Please doc me! */
  template<>
  class PrismQuadraturePoints<3>
  {
  public:
    enum { MAXP=6};
    enum { highest_order=2 };

    //! initialize quadrature points on the interval for all orders
    PrismQuadraturePoints ()
    {
      int m = 0;
      O[m] = 0;

      // polynom degree 0  ???
      m = 6;
      G[m][0][0] = 0.0;
      G[m][0][1] = 0.0;
      G[m][0][2] = 0.0;

      G[m][1][0] = 1.0;
      G[m][1][1] = 0.0;
      G[m][1][2] = 0.0;

      G[m][2][0] = 0.0;
      G[m][2][1] = 1.0;
      G[m][2][2] = 0.0;

      G[m][3][0] = 0.0;
      G[m][3][1] = 0.0;
      G[m][3][2] = 1.0;

      G[m][4][0] = 1.0;
      G[m][4][1] = 0.0;
      G[m][4][2] = 1.0;

      G[m][5][0] = 0.0;
      G[m][5][1] = 0.1;
      G[m][5][2] = 1.0;

      W[m][0] = 0.16666666666666666 / 2.0;
      W[m][1] = 0.16666666666666666 / 2.0;
      W[m][2] = 0.16666666666666666 / 2.0;
      W[m][3] = 0.16666666666666666 / 2.0;
      W[m][4] = 0.16666666666666666 / 2.0;
      W[m][5] = 0.16666666666666666 / 2.0;

      O[m] = 0;  // verify ????????


      // polynom degree 2  ???
      m = 6;
      G[m][0][0] =0.66666666666666666 ;
      G[m][0][1] =0.16666666666666666 ;
      G[m][0][2] =0.211324865405187 ;

      G[m][1][0] = 0.16666666666666666;
      G[m][1][1] =0.66666666666666666 ;
      G[m][1][2] = 0.211324865405187;

      G[m][2][0] = 0.16666666666666666;
      G[m][2][1] = 0.16666666666666666;
      G[m][2][2] = 0.211324865405187;

      G[m][3][0] = 0.66666666666666666;
      G[m][3][1] = 0.16666666666666666;
      G[m][3][2] = 0.788675134594813;

      G[m][4][0] = 0.16666666666666666;
      G[m][4][1] = 0.66666666666666666;
      G[m][4][2] = 0.788675134594813;

      G[m][5][0] = 0.16666666666666666;
      G[m][5][1] = 0.16666666666666666;
      G[m][5][2] = 0.788675134594813;

      W[m][0] = 0.16666666666666666 / 2.0;
      W[m][1] = 0.16666666666666666 / 2.0;
      W[m][2] = 0.16666666666666666 / 2.0;
      W[m][3] = 0.16666666666666666 / 2.0;
      W[m][4] = 0.16666666666666666 / 2.0;
      W[m][5] = 0.16666666666666666 / 2.0;

      O[m] = 2;  // verify ????????

    }

    /** \todo Please doc me! */
    FieldVector<double, 3> point(int m, int i)
    {
      return G[m][i];
    }

    /** \todo Please doc me! */
    double weight (int m, int i)
    {
      return W[m][i];
    }

    /** \todo Please doc me! */
    int order (int m)
    {
      return O[m];
    }

  private:
    FieldVector<double, 3> G[MAXP+1][MAXP]; //positions

    double W[MAXP+1][MAXP];     // weights associated with points
    int O[MAXP+1];              // order of the rule
  };


  /** \brief Singleton holding the Prism Quadrature points
     \ingroup Quadrature
   */
  template<int dim>
  struct PrismQuadraturePointsSingleton {
    static PrismQuadraturePoints<3> prqp;
  };

  /** \brief Singleton holding the Prism Quadrature points
     \ingroup Quadrature
   */
  template<>
  struct PrismQuadraturePointsSingleton<3> {
    static PrismQuadraturePoints<3> prqp;
  };

  /** \brief Quadrature rules for prisms
      \ingroup Quadrature
   */
  template<typename ct, int dim>
  class PrismQuadratureRule;

  /** \brief Quadrature rules for prisms
      \ingroup Quadrature
   */
  template<typename ct>
  class PrismQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {
  public:

    /** \brief The space dimension */
    enum { d = 3 };

    /** \brief The highest quadrature order available */
    enum { highest_order = 2 };

    ~PrismQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    PrismQuadratureRule(int p) : QuadratureRule<ct,3>(GeometryTypes::prism)
    {
      int m=6;
      this->delivered_order = PrismQuadraturePointsSingleton<3>::prqp.order(m);
      for(int i=0; i<m; ++i)
      {
        FieldVector<ct,3> local;
        for (int k=0; k<d; k++)
          local[k] = PrismQuadraturePointsSingleton<3>::prqp.point(m,i)[k];
        double weight =
          PrismQuadraturePointsSingleton<3>::prqp.weight(m,i);
        // put in container
        this->push_back(QuadraturePoint<ct,d>(local,weight));
      }
    }
  };

  /** \brief Factory class for creation of quadrature rules,
      depending on GeometryType, order and QuadratureType.

      The whole class is private and can only be accessed
      by the singleton container class QuadratureRules.
   */
  template<typename ctype, int dim>
  class QuadratureRuleFactory {
  private:
    friend class QuadratureRules<ctype, dim>;
    static unsigned maxOrder(const GeometryType &t, QuadratureType::Enum qt)
    {
      return TensorProductQuadratureRule<ctype,dim>::maxOrder(t.id(), qt);
    }
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      return TensorProductQuadratureRule<ctype,dim>(t.id(), p, qt);
    }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 0> {
  private:
    enum { dim = 0 };
    friend class QuadratureRules<ctype, dim>;
    static unsigned maxOrder(const GeometryType &t, QuadratureType::Enum qt)
    {
      if (t.isVertex())
      {
        return std::numeric_limits<int>::max();
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      if (t.isVertex())
      {
        return PointQuadratureRule<ctype>();
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 1> {
  private:
    enum { dim = 1 };
    friend class QuadratureRules<ctype, dim>;
    static unsigned maxOrder(const GeometryType &t, QuadratureType::Enum qt)
    {
      if (t.isLine())
      {
        switch (qt) {
        case QuadratureType::GaussLegendre :
          return GaussQuadratureRule1D<ctype>::highest_order;
        case QuadratureType::GaussJacobi_1_0 :
          return Jacobi1QuadratureRule1D<ctype>::highest_order;
        case QuadratureType::GaussJacobi_2_0 :
          return Jacobi2QuadratureRule1D<ctype>::highest_order;
        case QuadratureType::GaussLobatto :
          return GaussLobattoQuadratureRule1D<ctype>::highest_order;
        case QuadratureType::GaussJacobi_n_0 :
          return JacobiNQuadratureRule1D<ctype>::maxOrder();
        case QuadratureType::GaussRadauLeft :
          return GaussRadauLeftQuadratureRule1D<ctype>::highest_order;
        case QuadratureType::GaussRadauRight :
          return GaussRadauRightQuadratureRule1D<ctype>::highest_order;
        default :
          DUNE_THROW(Exception, "Unknown QuadratureType");
        }
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      if (t.isLine())
      {
        switch (qt) {
        case QuadratureType::GaussLegendre :
          return GaussQuadratureRule1D<ctype>(p);
        case QuadratureType::GaussJacobi_1_0 :
          return Jacobi1QuadratureRule1D<ctype>(p);
        case QuadratureType::GaussJacobi_2_0 :
          return Jacobi2QuadratureRule1D<ctype>(p);
        case QuadratureType::GaussLobatto :
          return GaussLobattoQuadratureRule1D<ctype>(p);
        case QuadratureType::GaussJacobi_n_0 :
          return JacobiNQuadratureRule1D<ctype>(p);
        case QuadratureType::GaussRadauLeft :
          return GaussRadauLeftQuadratureRule1D<ctype>(p);
        case QuadratureType::GaussRadauRight :
          return GaussRadauRightQuadratureRule1D<ctype>(p);
        default :
          DUNE_THROW(Exception, "Unknown QuadratureType");
        }
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 2> {
  private:
    enum { dim = 2 };
    friend class QuadratureRules<ctype, dim>;
    static unsigned maxOrder(const GeometryType &t, QuadratureType::Enum qt)
    {
      unsigned order =
        TensorProductQuadratureRule<ctype,dim>::maxOrder(t.id(), qt);
      if (t.isSimplex())
        order = std::max
          (order, unsigned(SimplexQuadratureRule<ctype,dim>::highest_order));
      return order;
    }
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {
      if (t.isSimplex()
        && ( qt == QuadratureType::GaussLegendre || qt == QuadratureType::GaussJacobi_n_0 )
        && p <= SimplexQuadratureRule<ctype,dim>::highest_order)
      {
        return SimplexQuadratureRule<ctype,dim>(p);
      }
      return TensorProductQuadratureRule<ctype,dim>(t.id(), p, qt);
    }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 3> {
  private:
    enum { dim = 3 };
    friend class QuadratureRules<ctype, dim>;
    static unsigned maxOrder(const GeometryType &t, QuadratureType::Enum qt)
    {
      unsigned order =
        TensorProductQuadratureRule<ctype,dim>::maxOrder(t.id(), qt);
      if (t.isSimplex())
        order = std::max
          (order, unsigned(SimplexQuadratureRule<ctype,dim>::highest_order));
      if (t.isPrism())
        order = std::max
          (order, unsigned(PrismQuadratureRule<ctype,dim>::highest_order));
      return order;
    }
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
    {

      if (t.isSimplex()
        && ( qt == QuadratureType::GaussLegendre || qt == QuadratureType::GaussJacobi_n_0 )
        && p <= SimplexQuadratureRule<ctype,dim>::highest_order)
      {
        return SimplexQuadratureRule<ctype,dim>(p);
      }
      if (t.isPrism()
        && qt == QuadratureType::GaussLegendre
        && p <= PrismQuadratureRule<ctype,dim>::highest_order)
      {
        return PrismQuadratureRule<ctype,dim>(p);
      }
      return TensorProductQuadratureRule<ctype,dim>(t.id(), p, qt);
    }
  };

} // end namespace

#endif // DUNE_GEOMETRY_QUADRATURERULES_HH
