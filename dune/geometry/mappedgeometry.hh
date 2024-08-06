// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_MAPPEDGEOMETRY_HH
#define DUNE_GEOMETRY_MAPPEDGEOMETRY_HH

#include <cassert>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>

#include <dune/common/copyableoptional.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/transpose.hh>
#include <dune/geometry/affinegeometry.hh> // for FieldMatrixHelper
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/utility/algorithms.hh>
#include <dune/geometry/utility/convergence.hh>

namespace Dune {

/**
 * \brief Geometry parametrized by a LocalFunction and a LocalGeometry.
 *
 * This class represents a geometry that is parametrized by the chained mapping
 * of a geometry `g` of type `Geo` and the given (differentiable) function `f`
 * of type `Map` as `f o g`.
 *
 * The geometry `g: LG -> GG` maps points from the local domain `LG` to its
 * global domain `GG`. It must fulfill the dune-grid geometry-concept. Examples
 * are the geometry of a grid element, the geometry of a \ref `ReferenceElement`
 * or any other geometry implementation in dune-geometry. The local domain `LG`
 * represents the local coordinates of the `MappedGeometry`.
 *
 * The function `f: LF -> GF` is a differentiable function in the sense of the
 * dune-functions concept. It maps coordinates from the global domain of the
 * geometry, `LF=GG` into the range type `GF` representing the global coordinates
 * of the `MappedGeometry`.
 *
 * \b Requirements:
 * Let `local` be of type `LG`, `g` of type `Geo` and `f` of type `Map`.
 * - `g` must be a model of `Concept::Geometry`
 * - `f` must be a model of `Concept::DifferentiableFunction<GF(LF)>`
 *
 * The following expressions must be valid:
 * - `GG gg = g.global(local)`: the "evaluation" of the geometry in its local domain.
 * - `GF gf = f(gg)`: the evaluation of the mapping `f` in its local domain that
 *      is the geometries' range.
 * - `auto df = derivative(f)`, `df(gg)`: The derivative of `f` w.r.t. its global
 *      coordinate `GF` evaluated at its local coordinated `LF=GG`.
 *
 * \tparam Map  Differentiable mapping `f: LF -> GF` with `LF = GG` the range of
 *              the geometry `g`.
 * \tparam Geo  A geometry type `g: LG -> GG` fulfilling the concept `Concept::Geometry`.
 **/
template <class Map, class Geo>
class MappedGeometry
{
public:
  /// type of local coordinates
  using LocalCoordinate = typename Geo::LocalCoordinate;

  /// type of global coordinates
  using GlobalCoordinate = std::remove_reference_t<decltype(std::declval<Map>()(std::declval<typename Geo::GlobalCoordinate>()))>;

  /// coordinate type
  using ctype = typename Geo::ctype;

  /// geometry dimension
  static constexpr int mydimension = LocalCoordinate::size();

  /// coordinate dimension
  static constexpr int coorddimension = GlobalCoordinate::size();

  /// type of volume
  using Volume = std::remove_reference_t<decltype(Dune::power(std::declval<ctype>(),mydimension))>;

  /// type of jacobian
  using Jacobian = FieldMatrix<ctype, coorddimension, mydimension>;

  /// type of jacobian transposed
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse
  using JacobianInverse = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

private:
  using ReferenceElements = Dune::ReferenceElements<ctype, mydimension>;
  using ReferenceElement = typename ReferenceElements::ReferenceElement;

protected:
  using MatrixHelper = Impl::FieldMatrixHelper<ctype>;

  // type of the mapping representation the geometry parametrization
  using Mapping = Map;

  // type of the geometry that is wrapped
  using Geometry = Geo;

  // type of a mapping representing the derivative of `Map` w.r.t. `GlobalCoordinate`
  using DerivativeMapping = std::remove_reference_t<decltype(derivative(std::declval<Map>()))>;

public:
  /**
   * \brief Constructor from mapping to parametrize the geometry.
   * \param[in]  mapping   A differentiable function for the parametrization of
   *                       the geometry
   * \param[in]  geometry  The geometry that is mapped
   * \param[in]  affine    Flag indicating whether the `MappedGeometry` represents
   *                       an affine mapping
   **/
  template <class Geo_, class Map_,
    std::enable_if_t<Dune::IsInteroperable<Map, Map_>::value, int> = 0,
    std::enable_if_t<Dune::IsInteroperable<Geo, Geo_>::value, int> = 0>
  MappedGeometry (Map_&& mapping, Geo_&& geometry, bool affine = false)
    : mapping_(std::forward<Map_>(mapping))
    , dMapping_(derivative(*mapping_))
    , geometry_(std::forward<Geo_>(geometry))
    , affine_(affine)
  {}

  /**
   * \brief Is this mapping affine? Not in general, since we don't know anything
   * about the mapping. The returned value can be given in the constructor of the
   * class.
   */
  bool affine () const
  {
    return affine_;
  }

  /// \brief Obtain the geometry type from the reference element.
  GeometryType type () const
  {
    return geometry_.type();
  }

  /// \brief Obtain number of corners of the corresponding reference element.
  int corners () const
  {
    return geometry_.corners();
  }

  /// \brief Obtain coordinates of the i-th corner.
  GlobalCoordinate corner (int i) const
  {
    assert( (i >= 0) && (i < corners()) );
    return mapping()(geometry_.corner(i));
  }

  /// \brief Map the center of the wrapped geometry.
  GlobalCoordinate center () const
  {
    return mapping()(geometry_.center());
  }

  /**
   *  \brief Evaluate the coordinate mapping.
   *
   *  Chained evaluation of the geometry and mapping from local coordinates in
   *  the reference element associated to the wrapped geometry.
   *
   *  \param[in] local  local coordinates
   *  \returns corresponding global coordinate
   **/
  GlobalCoordinate global (const LocalCoordinate& local) const
  {
    return mapping()(geometry_.global(local));
  }

  /**
   * \brief Evaluate the inverse coordinate mapping.
   * \param[in] y  Global coordinate to map
   * \param[in] opts  Parameters to control the behavior of the Gauss-Newton
   *                  algorithm.
   *
   * \return For given global coordinate `y` the local coordinate `x` that minimizes
   *         the function `(global(x) - y).two_norm2()` over the local coordinate
   *         space spanned by the reference element.
   *
   *  \throws In case of an error indicating that the local coordinate can not be
   *          obtained, a `Dune::Exception` is thrown, with an error code from
   *          \ref `GaussNewtonErrorCode`.
   *  \note It is not guaranteed that the resulting local coordinate is inside the
   *        reference element domain.
   **/
  LocalCoordinate local (const GlobalCoordinate& y, Impl::GaussNewtonOptions<ctype> opts = {}) const
  {
    LocalCoordinate x = refElement().position(0,0);
    Impl::GaussNewtonErrorCode err = Impl::gaussNewton(
      [&](const LocalCoordinate& local) { return this->global(local); },
      [&](const LocalCoordinate& local) { return this->jacobianTransposed(local); },
      y, x, opts
    );

    if (err != Impl::GaussNewtonErrorCode::OK)
      DUNE_THROW(Dune::Exception,
        "Local coordinate can not be recovered from global coordinate, error code = " << int(err) << "\n"
        << "  (global(x) - y).two_norm() = " << (global(x) - y).two_norm()
        << " > tol = " << opts.absTol);

    return x;
  }

  /**
   * \brief Obtain the integration element.
   *
   * If the Jacobian of the geometry is denoted by \f$J(x)\f$, the integration element
   * \f$\mu(x)\f$ is given by \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
   *
   * \param[in]  local  local coordinate to evaluate the integration element in.
   *
   * \returns the integration element \f$\mu(x)\f$.
   **/
  ctype integrationElement (const LocalCoordinate& local) const
  {
    return MatrixHelper::sqrtDetAAT(jacobianTransposed(local));
  }

  /**
   * \brief Obtain the volume of the mapping's image.
   *
   * Calculates the volume of the entity by numerical integration. Since the
   * polynomial order of the volume element is not known, iteratively compute
   * numerical integrals with increasing order of the quadrature rules, until
   * tolerance is reached.
   *
   * \param opts  An optional control over the convergence, providing a break
   *              tolerance and a maximal iteration count.
   **/
  Volume volume (Impl::ConvergenceOptions<ctype> opts = {}) const
  {
    Volume vol0 = volume(QuadratureRules<ctype, mydimension>::rule(type(), 1));
    if (affine())
      return vol0;

    using std::abs;
    for (int p = 2; p < opts.maxIt; ++p) {
      Volume vol1 = volume(QuadratureRules<ctype, mydimension>::rule(type(), p));
      if (abs(vol1 - vol0) < opts.absTol)
        return vol1;

      vol0 = vol1;
    }
    return vol0;
  }

  /// \brief Obtain the volume of the mapping's image by given quadrature rules.
  Volume volume (const QuadratureRule<ctype, mydimension>& quadRule) const
  {
    Volume vol(0);
    for (const auto& qp : quadRule)
      vol += integrationElement(qp.position()) * qp.weight();
    return vol;
  }

  /**
   * \brief Obtain the Jacobian.
   * \param[in]  local  local coordinate to evaluate Jacobian in
   * \returns the matrix corresponding to the Jacobian
   **/
  Jacobian jacobian (const LocalCoordinate& local) const
  {
    auto&& jLocal = geometry_.jacobian(local);
    auto&& jMapping = (*dMapping_)(geometry_.global(local));
    return jMapping * jLocal;
  }

  /**
   * \brief Obtain the transposed of the Jacobian.
   * \param[in]  local  local coordinate to evaluate Jacobian in
   * \returns the matrix corresponding to the transposed of the Jacobian
   **/
  JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    return transpose(jacobian(local));
  }

  /**
   * \brief Obtain the Jacobian's inverse.
   *
   * The Jacobian's inverse is defined as a pseudo-inverse. If we denote the
   * Jacobian by \f$J(x)\f$, the following condition holds:
   *   \f[ J^{-1}(x) J(x) = I. \f]
   **/
  JacobianInverse jacobianInverse (const LocalCoordinate& local) const
  {
    JacobianInverse out;
    MatrixHelper::leftInvA(jacobian(local), out);
    return out;
  }

  /**
   * \brief Obtain the transposed of the Jacobian's inverse.
   *
   * The Jacobian's inverse is defined as a pseudo-inverse. If we denote the
   * Jacobian by \f$J(x)\f$, the following condition holds:
   *   \f[ J^{-1}(x) J(x) = I. \f]
   **/
  JacobianInverseTransposed jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    return transpose(jacobianInverse(local));
  }

  /// \brief Obtain the reference-element related to this geometry.
  friend ReferenceElement referenceElement (const MappedGeometry& geometry)
  {
    return geometry.refElement();
  }

protected:
  // the internal stored reference element
  ReferenceElement refElement () const
  {
    return referenceElement(geometry_);
  }

private:
  // internal reference to the stored mapping
  const Mapping& mapping () const
  {
    return *mapping_;
  }

  // internal reference to the wrapped geometry
  const Geometry& geometry () const
  {
    return geometry_;
  }

private:
  /// Parametrization of the element
  CopyableOptional<Mapping> mapping_;

  /// The derivative of the \ref mapping_
  CopyableOptional<DerivativeMapping> dMapping_;

  /// The geometry that is wrapped
  Geometry geometry_;

  /// Flag whether the geometry mapping is affine
  bool affine_;
};

// deduction guides
template <class Map, class Geo>
MappedGeometry (const Map&, const Geo&)
  -> MappedGeometry<Map,Geo>;

template <class Map, class Geo>
MappedGeometry (const Map&, const Geo&, bool)
  -> MappedGeometry<Map,Geo>;

} // end namespace Dune

#endif // DUNE_GEOMETRY_MAPPEDGEOMETRY_HH
