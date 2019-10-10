#include <dune/geometry/quadraturerules.hh>
#include <functional>
#include <numeric>
#include <vector>

#include "quadraturerules.cc"

template<typename ct, int dim>
ct f(Dune::FieldVector<ct,dim> const& vec)
{
  auto const& x = vec[0];
  return x*x*x*x*x*x*x*x;
}

template<typename ct, int dim>
ct integrate(std::function<ct(Dune::FieldVector<ct,dim> const&)> const& f, Dune::QuadratureRule<ct,dim> const& rule)
{
  ct sum =0;
  auto applyFunctionAndSumUp = [f](ct const& sum, auto const& quadpoint) {
                       return sum + f(quadpoint.position())*quadpoint.weight();
                   };
  return std::accumulate(rule.begin(), rule.end(), sum, applyFunctionAndSumUp);
}
template<typename ct,int dim>
void printRule(Dune::QuadratureRule<ct,dim> const& rule)
{
  for(size_t i=0;i<rule.size();++i)
  {
    std::cout << "position = ";
    for(size_t j=0;j<rule[i].position().size();++j)
    {
      std::cout << rule[i].position()[j] << ",";
    }
    std::cout << "     weight = "<< rule[i].weight() << "\n";
  }
}
