// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_COMPOSITE_QUADRATURE_RULE_HH
#define DUNE_GEOMETRY_COMPOSITE_QUADRATURE_RULE_HH

/** \file
 * \brief Construct composite quadrature rules from other quadrature rules
 */

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/virtualrefinement.hh>

namespace Dune {

  /** \brief Construct composite quadrature rules from other quadrature rules
   *
   * \tparam ctype Type used for coordinates and quadrature weights
   * \tparam dim Dimension of the reference element
   */
  template <class ctype, int dim>
  class CompositeQuadratureRule
      : public Dune::QuadratureRule<ctype,dim>
  {
    public:
    /** \brief Construct composite quadrature rule
     * \param quad Base quadrature rule.  Element type of this rule must be simplex
     * \param refinement Number of uniform refinement steps
     */
    CompositeQuadratureRule(const Dune::QuadratureRule<ctype,dim>& quad, int refinement)
      : QuadratureRule<ctype,dim>(quad.type(), quad.order())
    {
      // Currently only works for simplices, because we are using the StaticRefinement
      assert(quad.type().isSimplex());

      typedef Dune::StaticRefinement<Dune::Impl::SimplexTopology<dim>::type::id,
                                     ctype,
                                     Dune::Impl::SimplexTopology<dim>::type::id,
                                     dim> Refinement;
      typedef typename Refinement::ElementIterator eIterator;

      ctype volume = Dune::ReferenceElements<ctype,dim>::general(quad.type()).volume();

      eIterator eSubEnd = Refinement::eEnd(refinement);
      eIterator eSubIt  = Refinement::eBegin(refinement);

      for (; eSubIt != eSubEnd; ++eSubIt) {

        // Percentage of the overall volume of this subelement
        ctype volumeFraction = eSubIt.geometry().volume() / volume;

        for (size_t i=0; i<quad.size(); i++) {

          this->push_back(Dune::QuadraturePoint<ctype,dim>(eSubIt.geometry().global(quad[i].position()),
                                                           volumeFraction*quad[i].weight()));

        }

      }

    }

  };

}

#endif   // DUNE_GEOMETRY_COMPOSITE_QUADRATURE_RULE_HH
