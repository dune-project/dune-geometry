#ifndef CALCULATE_SIMPLEX_RULES_H
#define CALCULATE_SIMPLEX_RULES_H

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
   * This class calculates quadrature rules of arbitrary order, that respect the
   * weightfunctions resulting from the conical product.
   * For more information see the comment in `tensorproductquadrature.hh`
   *
   */

  template<typename ct>
  class JacobiNQuadratureRule1D : public QuadratureRule<ct,1>
  {
  public:
    // compile time parameters
    enum { dim=1 };

    ~JacobiNQuadratureRule1D(){}
  //private:

    typedef QuadratureRule<ct, dim> Rule;

    friend class QuadratureRuleFactory<ct,dim>;

    template< class ctype, int dimension>
    friend class TensorProductQuadratureRule;

    JacobiNQuadratureRule1D (int const order, int const alpha=0)
      : Rule( GeometryTypes::line )
    {
      if (unsigned(order) > maxOrder())
        DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << order << " not supported!");
      auto&& rule = decideRule(order,alpha);
      for(unsigned int i=0; i< rule.size(); i++)
        this->push_back(rule[i]);
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

    template<typename type, std::enable_if_t<!(std::is_floating_point<type>::value)>* = nullptr>
    QuadratureRule<ct,1> UseLapackOrError( int const degree, int const alpha)
    {
      // TODO make an Error
      DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << degree << " not supported!");
    }

    template<typename type, std::enable_if_t<std::is_floating_point<type>::value>* = nullptr>
    QuadratureRule<ct,1> UseLapackOrError( int const degree, int const alpha)
    {
      return jacobiNodesWeights(degree, alpha);
    }

    QuadratureRule<ct,1> decideRuleExponent(int const degree, int const alpha)
    {
      switch(alpha)
      {
        case 0 :  return QuadratureRules<ct,1>::rule(GeometryTypes::line, degree, QuadratureType::GaussLegendre);
                  // the weights are not correct for the use of quadrature rules
        case 1 :  {
                    auto&& rule = QuadratureRules<ct,1>::rule(GeometryTypes::line, degree, QuadratureType::GaussJacobi_1_0);
                    QuadratureRule<ct,1> quadratureRule;
                    quadratureRule.reserve(rule.size());
                    for( auto qpoint : rule )
                      quadratureRule.push_back(QuadraturePoint<ct,1>(qpoint.position(),ct(0.5)*qpoint.weight()));
                    return quadratureRule;
                  }
        case 2 :  {
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

    // computes the nodes and weights for the weight function (1-x)^alpha, which are exact for given polynomials with degree "degree"
    QuadratureRule<ct,1> jacobiNodesWeights(int const degree, int const alpha)
    {
      using std::sqrt;

      typedef DynamicVector<ct> Vector;
      typedef DynamicMatrix<ct> Matrix;

      // compute the degree of the needed jacobi polymomial
      const int n = degree/2 +1;

      Matrix J(n,n,0);

      J[0][0] = (ct) -alpha/(2 + alpha);
      for(int i=1; i<n; ++i)
      {
        ct a_i = (ct) (2*i*(i + alpha)) /( (2*i + alpha - 1)*(2*i + alpha) );
        ct c_i = (ct) (2*i*(i + alpha)) /( (2*i + alpha + 1)*(2*i + alpha) );

        J[i][i] =  (ct) - alpha*alpha /( (2*i + alpha + 2)*(2*i + alpha) );
        J[i][i-1] = sqrt(a_i*c_i);
        J[i-1][i] = J[i][i-1];
      }

      DynamicVector<std::complex<double> > eigenValues(n,0);
      std::vector<DynamicVector<double> > eigenVectors(n, DynamicVector<double>(n,0));

      DynamicMatrixHelp::eigenValuesNonSym(J, eigenValues, &eigenVectors);

      ct mu = (ct) 1/(alpha + 1);

      QuadratureRule<ct,1> quadratureRule;
      quadratureRule.reserve(n);
      for (int i=0; i<n; ++i)
      {
        auto&& eV0 = eigenVectors[i][0];
        ct weight =  mu * eV0*eV0;
        Vector node(1,0.5*eigenValues[i].real() + 0.5);

        // bundle the nodes and the weights
        QuadraturePoint<ct,1> temp(node, weight);
        quadratureRule.push_back(temp);
      }

      return quadratureRule;

    }

  };

}

#endif
