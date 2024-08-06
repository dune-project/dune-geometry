// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_LOCALFINITEELEMENTGEOMETRY_HH
#define DUNE_GEOMETRY_LOCALFINITEELEMENTGEOMETRY_HH

#include <cassert>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/std/type_traits.hh>

#include <dune/geometry/affinegeometry.hh> // for FieldMatrixHelper
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/utility/algorithms.hh>
#include <dune/geometry/utility/convergence.hh>

namespace Dune {

/**
 * \brief Geometry implementation based on local-basis function parametrization.
 *
 * Parametrization of the geometry by any localfunction interpolated into a local
 * finite-element space.
 *
 * \tparam  LFE   Type of a local finite-element.
 * \tparam  cdim  Coordinate dimension.
 */
template <class LFE, int cdim>
class LocalFiniteElementGeometry
{
  using LocalFiniteElement = LFE;
  using LocalBasis = typename LFE::Traits::LocalBasisType;
  using LocalBasisTraits = typename LocalBasis::Traits;

public:
  /// coordinate type
  using ctype = typename LocalBasisTraits::DomainFieldType;

  /// geometry dimension
  static const int mydimension = LocalBasisTraits::dimDomain;

  /// coordinate dimension
  static const int coorddimension = cdim;

  /// type of local coordinates
  using LocalCoordinate = FieldVector<ctype, mydimension>;

  /// type of global coordinates
  using GlobalCoordinate = FieldVector<ctype, coorddimension>;

  /// type of volume
  using Volume = decltype(power(std::declval<ctype>(),mydimension));

  /// type of jacobian
  using Jacobian = FieldMatrix<ctype, coorddimension, mydimension>;

  /// type of jacobian transposed
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse
  using JacobianInverse = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

public:
  /// type of reference element
  using ReferenceElements = Dune::ReferenceElements<ctype, mydimension>;
  using ReferenceElement = typename ReferenceElements::ReferenceElement;

protected:
  using MatrixHelper = Impl::FieldMatrixHelper<ctype>;

public:
  /// \brief Default constructed geometry results in an empty/invalid representation.
  LocalFiniteElementGeometry () = default;

  /**
   * \brief Constructor from a vector of coefficients of the LocalBasis parametrizing
   *        the geometry.
   *
   * \param[in]  refElement  reference element for the geometry
   * \param[in]  localFE     Local finite-element to use for the parametrization
   * \param[in]  vertices    Coefficients of the local interpolation into the basis
   *
   * The vertices are the coefficients of a local interpolation in the local
   * finite-element. For Lagrange local bases, these correspond to vertices on
   * the curved geometry in the local Lagrange nodes.
   *
   * \note The vertices are stored internally, so if possible move an external
   *       vertex storage to this constructor
   **/
  LocalFiniteElementGeometry (const ReferenceElement& refElement,
                              const LocalFiniteElement& localFE,
                              std::vector<GlobalCoordinate> vertices)
    : refElement_(refElement)
    , localFE_(localFE)
    , vertices_(std::move(vertices))
  {
    assert(localFE_.size() == vertices_.size());
  }

  /**
   * \brief Constructor from a local parametrization function, mapping local to
   *        (curved) global coordinates.
   *
   * \param[in]  refElement      reference element for the geometry
   * \param[in]  localFE         Local finite-element to use for the parametrization
   * \param[in]  parametrization parametrization function with signature
   *                             `GlobalCoordinate(LocalCoordinate)`
   *
   * The parametrization function is not stored in the class, but interpolated into
   * the local finite-element basis and the computed interpolation coefficients
   * are stored.
   **/
  template <class Param,
    std::enable_if_t<std::is_invocable_r_v<GlobalCoordinate,Param,LocalCoordinate>, int> = 0>
  LocalFiniteElementGeometry (const ReferenceElement& refElement,
                              const LocalFiniteElement& localFE,
                              Param&& parametrization)
    : refElement_(refElement)
    , localFE_(localFE)
  {
    localFE_.localInterpolation().interpolate(parametrization, vertices_);
  }

  /**
   * \brief Constructor, forwarding to the other constructors that take a reference-element.
   *
   * \param[in]  gt       geometry type
   * \param[in]  args...  arguments passed to the other constructors
   **/
  template <class... Args>
  explicit LocalFiniteElementGeometry (GeometryType gt, Args&&... args)
    : LocalFiniteElementGeometry(ReferenceElements::general(gt), std::forward<Args>(args)...)
  {}

  /// \brief Obtain the (highest) polynomial order of the parametrization.
  int order () const
  {
    return localBasis().order();
  }

  /**
   * \brief Is this mapping affine?
   * Geometries of order 1 might be affine, but it needs to be checked on
   * non-simplex geometries.
   **/
  bool affine () const
  {
    if (!affine_)
      affine_.emplace(affineImpl());
    return *affine_;
  }

  /// \brief Obtain the name of the reference element.
  GeometryType type () const
  {
    return refElement_.type();
  }

  /// \brief Obtain number of corners of the corresponding reference element.
  int corners () const
  {
    return refElement_.size(mydimension);
  }

  /// \brief Obtain coordinates of the i-th corner.
  GlobalCoordinate corner (int i) const
  {
    assert( (i >= 0) && (i < corners()) );
    return global(refElement_.position(i, mydimension));
  }

  /// \brief Obtain the centroid of the mapping's image.
  GlobalCoordinate center () const
  {
    return global(refElement_.position(0, 0));
  }

  /**
   * \brief Evaluate the coordinate mapping.
   *
   * Implements a linear combination of local basis functions scaled by
   * the vertices as coefficients.
   *
   * \f[ global = \sum_i v_i \psi_i(local) \f]
   *
   * \param[in] local  local coordinate to map
   * \returns          corresponding global coordinate
   **/
  GlobalCoordinate global (const LocalCoordinate& local) const
  {
    thread_local std::vector<typename LocalBasisTraits::RangeType> shapeValues;
    localBasis().evaluateFunction(local, shapeValues);
    assert(shapeValues.size() == vertices_.size());

    GlobalCoordinate out(0);
    for (std::size_t i = 0; i < shapeValues.size(); ++i)
      out.axpy(shapeValues[i], vertices_[i]);

    return out;
  }

  /**
   * \brief Evaluate the inverse coordinate mapping.
   *
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
    LocalCoordinate x = refElement_.position(0,0);
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
    thread_local std::vector<typename LocalBasisTraits::JacobianType> shapeJacobians;
    localBasis().evaluateJacobian(local, shapeJacobians);
    assert(shapeJacobians.size() == vertices_.size());

    Jacobian out(0);
    for (std::size_t i = 0; i < shapeJacobians.size(); ++i) {
      for (int j = 0; j < Jacobian::rows; ++j) {
        shapeJacobians[i].umtv(vertices_[i][j], out[j]);
      }
    }
    return out;
  }

  /**
   * \brief Obtain the transposed of the Jacobian.
   * \param[in]  local  local coordinate to evaluate Jacobian in
   * \returns the matrix corresponding to the transposed of the Jacobian
   **/
  JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    return jacobian(local).transposed();
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
    return jacobianInverse(local).transposed();
  }

  /// \brief Obtain the reference-element related to this geometry.
  friend ReferenceElement referenceElement (const LocalFiniteElementGeometry& geometry)
  {
    return geometry.refElement_;
  }

  /// \brief Obtain the local finite-element.
  const LocalFiniteElement& finiteElement () const
  {
    return localFE_;
  }

  /// \brief Obtain the coefficients of the parametrization.
  const std::vector<GlobalCoordinate>& coefficients () const
  {
    return vertices_;
  }

  /// \brief The local basis of the stored local finite-element.
  const LocalBasis& localBasis () const
  {
    return localFE_.localBasis();
  }

private:

  bool affineImpl () const
  {
    if constexpr(mydimension == 0)
      // point geometries are always affine mappings
      return true;
    else {
      if (order() > 1)
        // higher-order parametrizations are by definition not affine
        return false;
      if constexpr(mydimension == 1)
        // linear line geometries are affine mappings
        return true;
      else {
        if (type().isSimplex())
          // linear simplex geometries are affine mappings
          return true;

        // multi-linear mappings on non-simplex geometries might be affine
        // as well. This is tested explicitly for all vertices by constructing
        // an affine mapping from dim+1 affine-independent corners and evaluating
        // at the other corners.
        auto refSimplex = referenceElement<ctype,mydimension>(GeometryTypes::simplex(mydimension));
        auto simplexIndices = Dune::range(refSimplex.size(mydimension));
        auto simplexCorners = Dune::transformedRangeView(simplexIndices,
          [&](int i) { return this->global(refSimplex.position(i,mydimension)); });
        AffineGeometry<ctype,mydimension,coorddimension> affineGeo(refSimplex,simplexCorners);
        using std::sqrt;
        const ctype tol = sqrt(std::numeric_limits<ctype>::epsilon());
        for (int i = 0; i < corners(); ++i) {
          const auto corner = refElement_.position(i,mydimension);
          if ((affineGeo.global(corner) - global(corner)).two_norm() > tol)
            return false;
        }
        return true;
      }
    }
  }

private:
  /// Reference of the geometry
  ReferenceElement refElement_{};

  /// A local finite-element
  LocalFiniteElement localFE_{};

  /// The (local finite-element) coefficients of the interpolating geometry
  std::vector<GlobalCoordinate> vertices_{};

  mutable std::optional<bool> affine_ = std::nullopt;
};

namespace Impl {

// extract the LocalCoordinate type from a LocalFiniteElement
template <class LFE>
using LocalCoordinate_t
  = FieldVector<typename LFE::Traits::LocalBasisType::Traits::DomainFieldType,
                LFE::Traits::LocalBasisType::Traits::dimDomain>;

} // end namespace Impl


// deduction guides
template <class I, class LFE, class GlobalCoordinate>
LocalFiniteElementGeometry (Geo::ReferenceElement<I>, const LFE&, std::vector<GlobalCoordinate>)
  -> LocalFiniteElementGeometry<LFE, GlobalCoordinate::dimension>;

template <class I, class LFE, class F,
          class Range = std::invoke_result_t<F,Impl::LocalCoordinate_t<LFE>>>
LocalFiniteElementGeometry (Geo::ReferenceElement<I>, const LFE&, const F&)
  -> LocalFiniteElementGeometry<LFE, Range::dimension>;

template <class LFE, class GlobalCoordinate>
LocalFiniteElementGeometry (GeometryType, const LFE& localFE, std::vector<GlobalCoordinate>)
  -> LocalFiniteElementGeometry<LFE, GlobalCoordinate::dimension>;

template <class LFE, class F,
          class Range = std::invoke_result_t<F,Impl::LocalCoordinate_t<LFE>>>
LocalFiniteElementGeometry (GeometryType, const LFE&, const F&)
  -> LocalFiniteElementGeometry<LFE, Range::dimension>;

} // namespace Dune

#endif // DUNE_GEOMETRY_LOCALFINITEELEMENTGEOMETRY_HH
