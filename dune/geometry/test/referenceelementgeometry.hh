// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TEST_REFERENCEELEMENTGEOMETRY_HH
#define DUNE_GEOMETRY_TEST_REFERENCEELEMENTGEOMETRY_HH

#include <type_traits>

#include <dune/geometry/referenceelements.hh>

namespace Dune {
namespace Impl {

// Placeholder type for a trivial identity matrix without any functionality
struct IdentityMatrix
{
  // multiply Id * A
  template <class A>
  friend const A& operator* (IdentityMatrix, const A& a) { return a; }

  // multiply A * Id
  template <class A>
  friend const A& operator* (const A& a, IdentityMatrix) { return a; }

  // multiply Id * Id
  friend IdentityMatrix operator* (IdentityMatrix, IdentityMatrix) { return {}; }

  friend std::ostream& operator<< (std::ostream& out, IdentityMatrix)
  {
    return out << "I";
  }

  // cast into FieldMatrix
  template <class K, int n>
  operator FieldMatrix<K,n,n> () const
  {
    FieldMatrix<K,n,n> I;
    for (int i = 0; i < n; ++i)
      I[i][i] = K(1);
    return I;
  }
};


/**
 * \brief Represent an identity map on a reference element as geometry.
 *
 * Each reference element `RefElem` provides an element-geometry that
 * is essentially a representation of the identity map. But it is typically
 * implemented as an \ref `AffineGeometry`. The `ReferenceElementGeometry`
 * instead uses exact identity with a special jacobian matrix encoded
 * in the `IdentityMatrix` class with optimized multiplications.
 *
 * The `ReferenceElementGeometry` inherits all functionality from
 * the geometry type of the reference element, except for `local()`,
 * `global()`, and `jacobinXXX()` methods.
 *
 * \tparam RefElem  The type of a reference element.
 **/
template <class RefElem>
class ReferenceElementGeometry
    : public RefElem::template Codim<0>::Geometry
{
  using Base = typename RefElem::template Codim<0>::Geometry;

public:
  using LocalCoordinate = typename RefElem::Coordinate;
  using GlobalCoordinate = typename RefElem::Coordinate;
  using Jacobian = Impl::IdentityMatrix;
  using JacobianTransposed = Impl::IdentityMatrix;
  using JacobianInverse = Impl::IdentityMatrix;
  using JacobianInverseTransposed = Impl::IdentityMatrix;

public:
  constexpr explicit ReferenceElementGeometry (const RefElem& refElem)
    : Base{refElem.template geometry<0>(0)}
  {}

  /// \brief Evaluate the inverse coordinate mapping, this is an identity mapping
  constexpr const LocalCoordinate& local (const GlobalCoordinate& global) const noexcept
  {
    return global;
  }

  /// \brief Evaluate the coordinate mapping, this is an identity mapping
  constexpr const GlobalCoordinate& global (const LocalCoordinate& local) const noexcept
  {
    return local;
  }

  /// \brief Obtain the Jacobian, this is an identity matrix.
  constexpr Jacobian jacobian (const LocalCoordinate& local) const noexcept
  {
    return Jacobian{};
  }

  /// \brief Obtain the transposed of the Jacobian, this is an identity matrix.
  constexpr JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const noexcept
  {
    return JacobianTransposed{};
  }

  /// \brief obtain the Jacobian's inverse, this is an identity matrix.
  constexpr JacobianInverse jacobianInverse (const LocalCoordinate& local) const noexcept
  {
    return JacobianInverse{};
  }

  /// \brief obtain the transposed of the Jacobian's inverse, this is an identity matrix.
  constexpr JacobianInverseTransposed jacobianInverseTransposed (const LocalCoordinate& local) const noexcept
  {
    return JacobianInverseTransposed{};
  }
};


/**
 * \brief Represent an identity map on a reference element with jacobians given by a real geometry.
 *
 * This geometry is a mixture of the `ReferenceElementGeometry` and the `Geometry`
 * given as a template parameter. Everything is taken from the `ReferenceElementGeometry`, that is,
 * essentially an identity mapping on the reference element, but the jacobian functions
 * are taken from the real `Geometry`.
 *
 * The application of the `LocalDerivativeGeometry` is inside a `MappendGeometry` together with a
 * local-function mapping from dune-functions that implements a derivative w.r.t. global coordinates,
 * instead of a derivative w.r.t. to local coordinates. A chaining of this function with a
 * `ReferenceElementGeometry` would need an additional transformation of the Jacobians that is
 * circumvented by using this `LocalDerivativeGeometry` directly.
 *
 * \tparam Geometry  The type of a real element geometry providing the `jacobianXXX` methods.
 **/
template <class Geometry>
class LocalDerivativeGeometry
    : public ReferenceElementGeometry<
        typename Dune::ReferenceElements<typename Geometry::ctype, Geometry::mydimension>::ReferenceElement >
{
  using ReferenceElements = Dune::ReferenceElements<typename Geometry::ctype, Geometry::mydimension>;
  using ReferenceElement = typename ReferenceElements::ReferenceElement;
  using Base = ReferenceElementGeometry<ReferenceElement>;

public:
  using LocalCoordinate = typename Geometry::LocalCoordinate;
  using Jacobian = typename Geometry::Jacobian;
  using JacobianTransposed = typename Geometry::JacobianTransposed;
  using JacobianInverse = typename Geometry::JacobianInverse;
  using JacobianInverseTransposed = typename Geometry::JacobianInverseTransposed;

public:
  explicit LocalDerivativeGeometry (const Geometry& geometry) noexcept
    : Base(referenceElement(geometry))
    , geometry_(geometry)
  {}

  /// \brief Obtain the Jacobian
  Jacobian jacobian (const LocalCoordinate& local) const
  {
    return geometry_.jacobian(local);
  }

  /// \brief Obtain the transposed of the Jacobian
  JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    return geometry_.jacobianTransposed(local);
  }

  /// \brief obtain the Jacobian's inverse
  JacobianInverse jacobianInverse (const LocalCoordinate& local) const
  {
    return geometry_.jacobianInverse(local);
  }

  /// \brief obtain the transposed of the Jacobian's inverse
  JacobianInverseTransposed jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    return geometry_.jacobianInverseTransposed(local);
  }

private:
  Geometry geometry_;
};

} // end namespace Impl
} // end namespace Dune

#endif // DUNE_GEOMETRY_TEST_REFERENCEELEMENTGEOMETRY_HH
