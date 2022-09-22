// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_QUADRATURERULES_JACOBI_N_0_H
#define DUNE_GEOMETRY_QUADRATURERULES_JACOBI_N_0_H

#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrixev.hh>


namespace Dune {

  /*
   * This class calculates 1D quadrature rules of dynamic order, that respect the
   * weightfunctions resulting from the conical product.
   * For quadrature rules with order <= highest order (61) and dimension <= 3 exact quadrature rules are used.
   * Else the quadrature rules are computed by Lapack and thus a floating point type is required.
   * For more information see the comment in `tensorproductquadrature.hh`
   *
   */

  template<typename ct, int dim>
  class JacobiNQuadratureRule;

  template<typename ct>
  using JacobiNQuadratureRule1D = JacobiNQuadratureRule<ct,1>;

  template<typename ct>
  class JacobiNQuadratureRule<ct,1> : public QuadratureRule<ct,1>
  {
  public:
    // compile time parameters
    constexpr static int dim = 1;

  private:

    typedef QuadratureRule<ct, dim> Rule;

    friend class QuadratureRuleFactory<ct,dim>;

    template< class ctype, int dimension>
    friend class TensorProductQuadratureRule;

    explicit JacobiNQuadratureRule (int const order, int const alpha=0)
      : Rule( GeometryTypes::line )
    {
      if (unsigned(order) > maxOrder())
        DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << order << " not supported!");
      auto&& rule = decideRule(order,alpha);
      for( auto qpoint : rule )
        this->push_back(qpoint);
      this->delivered_order = 2*rule.size()-1;

    }

    static unsigned maxOrder()
    {
      return 127; // can be changed
    }

    QuadratureRule<ct,1> decideRule(int const degree, int const alpha)
    {
      const auto maxOrder = std::min(   unsigned(GaussQuadratureRule1D<ct>::highest_order),
                              std::min( unsigned(Jacobi1QuadratureRule1D<ct>::highest_order),
                                        unsigned(Jacobi2QuadratureRule1D<ct>::highest_order))
                              );
      return unsigned(degree) < maxOrder ? decideRuleExponent(degree,alpha) : UseLapackOrError<ct>(degree, alpha);
    }

#if HAVE_LAPACK
    template<typename type, std::enable_if_t<!(std::is_floating_point<type>::value)>* = nullptr>
    QuadratureRule<ct,1> UseLapackOrError( int const degree, int const alpha)
    {
      DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule with degree: "<< degree << " and jacobi exponent: "<< alpha<< " is not supported for this type!");
    }

    template<typename type, std::enable_if_t<std::is_floating_point<type>::value>* = nullptr>
    QuadratureRule<ct,1> UseLapackOrError( int const degree, int const alpha)
    {
      return jacobiNodesWeights<ct>(degree, alpha);
    }
#else
    template<typename type>
    QuadratureRule<ct,1> UseLapackOrError( int const, int const)
    {
      DUNE_THROW(NotImplemented, "LAPACK must be enable to use JacobiN quadrature rules.");
    }
#endif

    QuadratureRule<ct,1> decideRuleExponent(int const degree, int const alpha)
    {
      switch(alpha)
      {
        case 0 :  return QuadratureRules<ct,1>::rule(GeometryTypes::line, degree, QuadratureType::GaussLegendre);
        case 1 :  {
                    // scale already existing weights by 0.5
                    auto&& rule = QuadratureRules<ct,1>::rule(GeometryTypes::line, degree, QuadratureType::GaussJacobi_1_0);
                    QuadratureRule<ct,1> quadratureRule;
                    quadratureRule.reserve(rule.size());
                    for( auto qpoint : rule )
                      quadratureRule.push_back(QuadraturePoint<ct,1>(qpoint.position(),ct(0.5)*qpoint.weight()));
                    return quadratureRule;
                  }
        case 2 :  {
                    // scale already existing weights by 0.25
                    auto&& rule = QuadratureRules<ct,1>::rule(GeometryTypes::line, degree, QuadratureType::GaussJacobi_2_0);
                    QuadratureRule<ct,1> quadratureRule;
                    quadratureRule.reserve(rule.size());
                    for( auto qpoint : rule )
                      quadratureRule.push_back(QuadraturePoint<ct,1>(qpoint.position(),ct(0.25)*qpoint.weight()));
                    return quadratureRule;
                  }
        default : return UseLapackOrError<ct>(degree,alpha);
      }
    }

#if HAVE_LAPACK
    // computes the nodes and weights for the weight function (1-x)^alpha, which are exact for given polynomials with degree "degree"
    template<typename ctype, std::enable_if_t<std::is_floating_point<ctype>::value>* = nullptr>
    QuadratureRule<ct,1> jacobiNodesWeights(int const degree, int const alpha)
    {
      using std::sqrt;

      // compute the degree of the needed jacobi polynomial
      const int n = degree/2 +1;

      DynamicMatrix<double> J(n,n,0);

      J[0][0] = -double(alpha)/(2 + alpha);
      for(int i=1; i<n; ++i)
      {
        ct a_i = double(2*i*(i + alpha)) /( (2*i + alpha - 1)*(2*i + alpha) );
        ct c_i = double(2*i*(i + alpha)) /( (2*i + alpha + 1)*(2*i + alpha) );

        J[i][i] =  -double(alpha*alpha) /( (2*i + alpha + 2)*(2*i + alpha) );
        J[i][i-1] = sqrt(a_i*c_i);
        J[i-1][i] = J[i][i-1];
      }

      DynamicVector<std::complex<double> > eigenValues(n,0);
      std::vector<DynamicVector<double> > eigenVectors(n, DynamicVector<double>(n,0));

      DynamicMatrixHelp::eigenValuesNonSym(J, eigenValues, &eigenVectors);

      double mu = 1.0/(alpha + 1);

      QuadratureRule<ct,1> quadratureRule;
      quadratureRule.reserve(n);
      for (int i=0; i<n; ++i)
      {
        auto&& eV0 = eigenVectors[i][0];
        ct weight =  mu * eV0*eV0;
        DynamicVector<ct> node(1,0.5*eigenValues[i].real() + 0.5);

        // bundle the nodes and the weights
        QuadraturePoint<ct,1> temp(node, weight);
        quadratureRule.push_back(temp);
      }

      return quadratureRule;

    }
#endif

  };

}

#endif
