// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TEST_LOCALFINITEELEMENT_HH
#define DUNE_GEOMETRY_TEST_LOCALFINITEELEMENT_HH

#include <algorithm>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

namespace Dune::Impl {

template <class D, class R, unsigned int dim>
struct ScalarLocalBasisTraits
{
  using DomainFieldType = D;
  using RangeFieldType = R;
  using DomainType = FieldVector<D,dim>;
  using RangeType = FieldVector<R,1>;
  using JacobianType = FieldMatrix<R,1,dim>;

  enum {
    dimDomain = dim,
    dimRange = 1
  };
};

/// Lagrange shape functions of order 1 on the reference simplex
template <class D, class R, unsigned int dim>
class P1LocalBasis
{
public:
  using Traits = ScalarLocalBasisTraits<D,R,dim>;

  /// Number of shape functions
  static constexpr unsigned int size () { return dim+1; }

  /// Evaluate all shape functions
  void evaluateFunction (const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(size());
    out[0] = 1.0;
    for (unsigned int i=0; i<dim; i++)
    {
      out[0]  -= x[i];
      out[i+1] = x[i];
    }
  }

  /// Evaluate Jacobian of all shape functions
  void evaluateJacobian (const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
  {
    out.resize(size());
    std::fill(out[0][0].begin(), out[0][0].end(), -1);

    for (unsigned int i=0; i<dim; i++)
      for (unsigned int j=0; j<dim; j++)
        out[i+1][0][j] = (i==j);
  }

  /// Polynomial order of the shape functions
  static constexpr unsigned int order () {  return 1; }
};

/// Lagrange shape functions of order 1 on the reference cube
template <class D, class R, unsigned int dim>
class Q1LocalBasis
{
public:
  using Traits = ScalarLocalBasisTraits<D,R,dim>;

  /// Number of shape functions
  static constexpr unsigned int size () { return power(2, dim); }

  /// Evaluate all shape functions
  void evaluateFunction (const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(size());
    for (size_t i=0; i<size(); i++)
    {
      out[i] = 1;

      for (unsigned int j=0; j<dim; j++)
        // if j-th bit of i is set multiply with x[j], else with 1-x[j]
        out[i] *= (i & (1<<j)) ? x[j] :  1-x[j];
    }
  }

  /// Evaluate Jacobian of all shape functions
  void evaluateJacobian (const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
  {
    out.resize(size());

    // Loop over all shape functions
    for (unsigned int i=0; i<size(); i++)
    {
      // Loop over all coordinate directions
      for (unsigned int j=0; j<dim; j++)
      {
        // Initialize: the overall expression is a product
        // if j-th bit of i is set to 1, else -11
        out[i][0][j] = (i & (1<<j)) ? 1 : -1;

        for (unsigned int l=0; l<dim; l++)
        {
          if (j!=l)
            // if l-th bit of i is set multiply with x[l], else with 1-x[l]
            out[i][0][j] *= (i & (1<<l)) ? x[l] :  1-x[l];
        }
      }
    }
  }

  /// Polynomial order of the shape functions
  static constexpr unsigned int order () { return 1; }
};


template <class LB>
class P1LocalInterpolation
{
public:
  /// Evaluate a given function at the Lagrange nodes
  template <class F, class C>
  void interpolate (F f, std::vector<C>& out) const
  {
    constexpr auto dim = LB::Traits::dimDomain;
    out.resize(LB::size());

    // vertex 0
    typename LB::Traits::DomainType x;
    std::fill(x.begin(), x.end(), 0);
    out[0] = f(x);

    // remaining vertices
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++)
        x[j] = (i==j);

      out[i+1] = f(x);
    }
  }
};

template <class LB>
class Q1LocalInterpolation
{
public:
  /// Evaluate a given function at the Lagrange nodes
  template <class F, class C>
  void interpolate (F f, std::vector<C>& out) const
  {
    constexpr auto dim = LB::Traits::dimDomain;
    out.resize(LB::size());

    typename LB::Traits::DomainType x;
    for (unsigned int i=0; i<LB::size(); i++)
    {
      // Generate coordinate of the i-th corner of the reference cube
      for (int j=0; j<dim; j++)
        x[j] = (i & (1<<j)) ? 1.0 : 0.0;

      out[i] = f(x);
    }
  }
};


/// Wrapper for local basis and local interpolation
template <class LB, template <class> class LI>
class LocalFiniteElement
{
public:
  struct Traits
  {
    using LocalBasisType = LB;
    using LocalInterpolationType = LI<LB>;
  };

  const auto& localBasis () const { return basis_; }
  const auto& localInterpolation () const { return interpolation_; }

  /// The number of shape functions
  static constexpr std::size_t size () { return LB::size(); }

private:
  LB basis_;
  LI<LB> interpolation_;
};

template <class D, class R, int d>
using P1LocalFiniteElement = LocalFiniteElement<P1LocalBasis<D,R,d>, P1LocalInterpolation>;

template <class D, class R, int d>
using Q1LocalFiniteElement = LocalFiniteElement<Q1LocalBasis<D,R,d>, Q1LocalInterpolation>;



template <class LFE, int cdim,
          class R = typename LFE::Traits::LocalBasisType::Traits::RangeFieldType>
class LocalFiniteElementFunction
{
  using LocalFiniteElement = LFE;
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using LocalBasisRange = typename LocalBasis::Traits::RangeType;
  using LocalBasisJacobian = typename LocalBasis::Traits::JacobianType;
  using Domain = typename LocalBasis::Traits::DomainType;
  using Range = FieldVector<R,cdim>;
  using Jacobian = FieldMatrix<R,cdim,LocalBasis::Traits::dimDomain>;

  static_assert(LocalBasis::Traits::dimRange == 1);

public:
  LocalFiniteElementFunction () = default;
  LocalFiniteElementFunction (const LocalFiniteElement& lfe, std::vector<Range> coefficients)
    : lfe_(lfe)
    , coefficients_(std::move(coefficients))
  {}

  Range operator() (const Domain& local) const
  {
    thread_local std::vector<LocalBasisRange> shapeValues;
    lfe_.localBasis().evaluateFunction(local, shapeValues);
    assert(shapeValues.size() == coefficients_.size());
    Range range(0);
    for (std::size_t i = 0; i < shapeValues.size(); ++i)
      range.axpy(shapeValues[i], coefficients_[i]);
    return range;
  }

  friend auto derivative (const LocalFiniteElementFunction& f)
  {
    return [&lfe=f.lfe_,coefficients=f.coefficients_](const Domain& local) -> Jacobian
    {
      thread_local std::vector<LocalBasisJacobian> shapeJacobians;
      lfe.localBasis().evaluateJacobian(local, shapeJacobians);
      assert(shapeJacobians.size() == coefficients.size());
      Jacobian jacobian(0);
      for (std::size_t i = 0; i < shapeJacobians.size(); ++i) {
        for (int j = 0; j < Jacobian::rows; ++j) {
          shapeJacobians[i].umtv(coefficients[i][j], jacobian[j]);
        }
      }
      return jacobian;
    };
  }

private:
  LocalFiniteElement lfe_{};
  std::vector<Range> coefficients_{};
};

} // end namespace Dune::Impl

#endif // DUNE_GEOMETRY_TEST_LOCALFINITEELEMENT_HH
