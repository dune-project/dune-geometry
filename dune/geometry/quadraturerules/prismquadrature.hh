// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_QUADRATURE_PRISM_HH
#define DUNE_GEOMETRY_QUADRATURE_PRISM_HH

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif

namespace Dune {

  /***********************************
   * quadrature for Prism
   **********************************/

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
    /** \brief The highest quadrature order available */
    constexpr static int highest_order = 2;

  private:
    friend class QuadratureRuleFactory<ct,3>;
    PrismQuadratureRule(int p);
    ~PrismQuadratureRule(){}
  };

  /** \todo Please doc me! */
  template<int dim>
  class PrismQuadraturePoints;

  /** \brief Singleton holding the Prism Quadrature points
     \ingroup Quadrature
   */
  template<int dim>
  struct PrismQuadraturePointsSingleton {};

  /** \brief Singleton holding the Prism Quadrature points
     \ingroup Quadrature
   */
  template<>
  struct PrismQuadraturePointsSingleton<3> {
    static PrismQuadraturePoints<3> prqp;
  };

  /** \todo Please doc me! */
  template<>
  class PrismQuadraturePoints<3>
  {
  public:
    constexpr static int MAXP = 6;
    constexpr static int highest_order = 2;

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

  template<typename ct>
  PrismQuadratureRule<ct,3>::PrismQuadratureRule(int /* p */) : QuadratureRule<ct,3>(GeometryTypes::prism)
  {
    int m=6;
    this->delivered_order = PrismQuadraturePointsSingleton<3>::prqp.order(m);
    for(int i=0; i<m; ++i)
    {
      FieldVector<ct,3> local = PrismQuadraturePointsSingleton<3>::prqp.point(m,i);
      ct weight = PrismQuadraturePointsSingleton<3>::prqp.weight(m,i);
      // put in container
      this->push_back(QuadraturePoint<ct,3>(local,weight));
    }
  }

} // namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_PRISM_HH
