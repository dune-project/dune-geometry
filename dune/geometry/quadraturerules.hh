// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

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
    constexpr static int dimension = dim;

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
      *                     In this case the 1D rules eventually need higher order to compensate occurring weight functions(i.e. simplices).
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
      *                     In this case the 1D rules respect eventually occurring weight functions(i.e. simplices).
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
    constexpr static int d = dim;

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
    using QuadratureRule = Dune::QuadratureRule<ctype, dim>;

    // indexed by quadrature order
    using QuadratureOrderVector = std::vector<std::pair<std::once_flag, QuadratureRule> >;

    // indexed by geometry type
    using GeometryTypeVector = std::vector<std::pair<std::once_flag, QuadratureOrderVector> >;

    // indexed by quadrature type enum
    using QuadratureCacheVector = std::vector<std::pair<std::once_flag, GeometryTypeVector> >;

    //! real rule creator
    DUNE_EXPORT const QuadratureRule& _rule(const GeometryType& t, int p, QuadratureType::Enum qt=QuadratureType::GaussLegendre) const
    {
      assert(t.dim()==dim);

      DUNE_ASSERT_CALL_ONCE();

      static QuadratureCacheVector quadratureCache(QuadratureType::size);

      auto& [ onceFlagQuadratureType, geometryTypes ] = quadratureCache[qt];
      // initialize geometry types for this quadrature type once
      std::call_once(onceFlagQuadratureType, [&types = geometryTypes]{
        types = GeometryTypeVector(LocalGeometryTypeIndex::size(dim));
      });

      auto& [ onceFlagGeometryType, quadratureOrders ] = geometryTypes[LocalGeometryTypeIndex::index(t)];
      // initialize quadrature orders for this geometry type and quadrature type once
      std::call_once(onceFlagGeometryType, [&, &orders = quadratureOrders]{
        // we only need one quadrature rule for points, not maxint
        const auto numRules = dim == 0 ? 1 : QuadratureRuleFactory<ctype,dim>::maxOrder(t, qt)+1;
        orders = QuadratureOrderVector(numRules);
      });

      // we only have one quadrature rule for points
      auto& [ onceFlagQuadratureOrder, quadratureRule ] = quadratureOrders[dim == 0 ? 0 : p];
      // initialize quadrature rule once
      std::call_once(onceFlagQuadratureOrder, [&, &rule = quadratureRule]{
        rule = QuadratureRuleFactory<ctype,dim>::rule(t, p, qt);
      });

      return quadratureRule;
    }

    //! singleton provider
    DUNE_EXPORT static QuadratureRules& instance()
    {
      static QuadratureRules instance;
      return instance;
    }

    //! private constructor
    QuadratureRules () = default;
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

#define DUNE_INCLUDING_IMPLEMENTATION

// 0d rules
#include "quadraturerules/pointquadrature.hh"
// 1d rules
#include "quadraturerules/gausslobattoquadrature.hh"
#include "quadraturerules/gaussquadrature.hh"
#include "quadraturerules/gaussradauleftquadrature.hh"
#include "quadraturerules/gaussradaurightquadrature.hh"
#include "quadraturerules/jacobi1quadrature.hh"
#include "quadraturerules/jacobi2quadrature.hh"
#include "quadraturerules/jacobiNquadrature.hh"
// 3d rules
#include "quadraturerules/prismquadrature.hh"
// general rules
#include "quadraturerules/simplexquadrature.hh"
#include "quadraturerules/tensorproductquadrature.hh"

#undef DUNE_INCLUDING_IMPLEMENTATION

namespace Dune {

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
    constexpr static int dim = 0;
    friend class QuadratureRules<ctype, dim>;
    static unsigned maxOrder(const GeometryType &t, QuadratureType::Enum)
    {
      if (t.isVertex())
      {
        return std::numeric_limits<int>::max();
      }
      DUNE_THROW(Exception, "Unknown GeometryType");
    }
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int , QuadratureType::Enum)
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
    constexpr static int dim = 1;
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
    constexpr static int dim = 2;
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
    constexpr static int dim = 3;
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

#ifndef DUNE_NO_EXTERN_QUADRATURERULES
  extern template class GaussLobattoQuadratureRule<double, 1>;
  extern template class GaussQuadratureRule<double, 1>;
  extern template class GaussRadauLeftQuadratureRule<double, 1>;
  extern template class GaussRadauRightQuadratureRule<double, 1>;
  extern template class Jacobi1QuadratureRule<double, 1>;
  extern template class Jacobi2QuadratureRule<double, 1>;
  extern template class JacobiNQuadratureRule<double, 1>;
  extern template class PrismQuadratureRule<double, 3>;
  extern template class SimplexQuadratureRule<double, 2>;
  extern template class SimplexQuadratureRule<double, 3>;
#endif // !DUNE_NO_EXTERN_QUADRATURERULES

} // end namespace

#endif // DUNE_GEOMETRY_QUADRATURERULES_HH
