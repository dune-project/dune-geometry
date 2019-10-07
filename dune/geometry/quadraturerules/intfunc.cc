#include "config.h"
#include <iostream>

#include "intfunc.hh"
#include "jacobi_n_0.hh"

#include <dune/geometry/quadraturerules.hh>

int main()
{
  using namespace Dune;
  const int dim = 5;
  const int order = 9;

  std::cout << "  dim = "<<dim <<", order = " <<order<<" \n";
  // get rule
  const auto& quadRuleOld = QuadratureRules<double,dim>::rule(GeometryType::simplex, order);

  std::cout << "--------------------\nquadRuleOld "<<quadRuleOld.size() <<" = \n";
  //printRule(quadRuleOld);
  std::cout << "int(f) = "<< integrate<double,dim>(f<double,dim>, quadRuleOld) <<"\n";

  // get new rule
  const auto& quadRule = QuadratureRules<double,dim>::rule(GeometryType::simplex, order, Dune::QuadratureType::GaussJacobi_n_0);
  std::cout << "--------------------\nQuadRule "<<quadRule.size() <<" = \n";
  //printRule(quadRule);

  std::cout << "int(f) = "<< integrate<double,dim>(f<double,dim>, quadRule) <<"\n";

  std::cout << "\n----------------------------------------------------------\n";
  const int alpha =1;
  std::cout << "jacobi alpha ="<< alpha <<"\n";
  const auto& rule = JacobiNQuadratureRule1D<double>(order, alpha);
  const auto& rule1 = JacobiNQuadratureRule1D<int>(order, alpha);
}
