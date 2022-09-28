// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_GEOMETRY_AXISALIGNED_CUBE_GEOMETRY_HH
#define DUNE_GEOMETRY_AXISALIGNED_CUBE_GEOMETRY_HH

/** \file
    \brief A geometry implementation for axis-aligned hypercubes
 */

#include <bitset>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/diagonalmatrix.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>


namespace Dune {

  /** \brief A geometry implementation for axis-aligned hypercubes
   *
   * This code is much faster than a generic implementation for hexahedral elements.
   * All methods use the fact that a geometry for axis-aligned cubes is basically just
   * a(n affine) scaling in the coordinate directions.
   *
   * If dim < coorddim then local coordinates need to be suitably mapped to global ones.
   * AxisAlignedCubeGeometry uses a special std::bitset 'axes' for this.  'axes' has coorddim
   * entries, of which precisely 'dim' need to be set.  Each set entry marks a local
   * coordinate, i.e., a coordinate in which the cube has extension.  The cube is flat
   * in all other directions.  Its coordinates in these directions is taking from
   * the array called 'lower', which specifies the lower left corner of the hypercube.
   *
   * In the case of dim==coorddim, the code goes into overdrive.  Then special code path's
   * are taken (statically) which omit the conditionals needed to sort out the embedding
   * of local into global coordinates.  Aggressive compiler/scheduler optimization becomes
   * possible.  Additionally, the types returned by the methods jacobianTransposed
   * and jacobianInverseTransposed are dedicated types for diagonal matrices (DiagonalMatrix).
   *
   * \tparam CoordType Type used for single coordinate coefficients
   * \tparam dim Dimension of the cube
   * \tparam coorddim Dimension of the space that the cube lives in
   */
  template <class CoordType, unsigned int dim, unsigned int coorddim>
  class AxisAlignedCubeGeometry
  {


  public:

    /** \brief Dimension of the cube element */
    constexpr static int mydimension = dim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    constexpr static int coorddimension = coorddim;

    /** \brief Type used for single coordinate coefficients */
    typedef CoordType ctype;

    /** \brief Type used for a vector of element coordinates */
    typedef FieldVector<ctype,dim> LocalCoordinate;

    /** \brief Type used for a vector of world coordinates */
    typedef FieldVector<ctype,coorddim> GlobalCoordinate;

    /** \brief Type used for volume */
    typedef ctype Volume;

    /** \brief Return type of jacobianTransposed

        This is a fast DiagonalMatrix if dim==coorddim, and a FieldMatrix otherwise.
        The FieldMatrix will never contain more than one entry per row,
        hence it could be replaced by something more efficient.
     */
    typedef typename std::conditional<dim==coorddim,
        DiagonalMatrix<ctype,dim>,
        FieldMatrix<ctype,dim,coorddim> >::type JacobianTransposed;

    /** \brief Return type of jacobianInverseTransposed

        This is a fast DiagonalMatrix if dim==coorddim, and a FieldMatrix otherwise.
        The FieldMatrix will never contain more than one entry per column,
        hence it could be replaced by something more efficient.
     */
    typedef typename std::conditional<dim==coorddim,
        DiagonalMatrix<ctype,dim>,
        FieldMatrix<ctype,coorddim,dim> >::type JacobianInverseTransposed;

    /**
     * \brief Return type of jacobian
     *
     * This is a fast DiagonalMatrix if dim==coorddim, and a FieldMatrix otherwise.
     * The FieldMatrix will never contain more than one entry per row,
     * hence it could be replaced by something more efficient.
     */
    using Jacobian = std::conditional_t<dim==coorddim, DiagonalMatrix<ctype,dim>, FieldMatrix<ctype,coorddim,dim> >;

    /**
     * \brief Return type of jacobianInverse
     *
     * This is a fast DiagonalMatrix if dim==coorddim, and a FieldMatrix otherwise.
     * The FieldMatrix will never contain more than one entry per row,
     * hence it could be replaced by something more efficient.
     */
    using JacobianInverse = std::conditional_t<dim==coorddim, DiagonalMatrix<ctype,dim>, FieldMatrix<ctype,dim,coorddim> >;

    /** \brief Constructor from a lower left and an upper right corner

        \note Only for dim==coorddim
     */
    AxisAlignedCubeGeometry(const Dune::FieldVector<ctype,coorddim> lower,
                            const Dune::FieldVector<ctype,coorddim> upper)
      : lower_(lower),
        upper_(upper),
        axes_()
    {
      static_assert(dim==coorddim, "Use this constructor only if dim==coorddim!");
      // all 'true', but is never actually used
      axes_ = (1<<coorddim)-1;
    }

    /** \brief Constructor from a lower left and an upper right corner
     *
     *  \param lower Coordinates for the lower left corner.
     *  \param upper Coordinates for the upper right corner.
     *  \param axes Each bit set to 'true' here corresponds to a local coordinate axis.
     *         In other words, precisely 'dim' bits must be set here.
     */
    AxisAlignedCubeGeometry(const Dune::FieldVector<ctype,coorddim> lower,
                            const Dune::FieldVector<ctype,coorddim> upper,
                            const std::bitset<coorddim>& axes)
      : lower_(lower),
        upper_(upper),
        axes_(axes)
    {
      assert(axes.count()==dim);
      for (size_t i=0; i<coorddim; i++)
        if (not axes_[i])
          upper_[i] = lower_[i];
    }

    /** \brief Constructor from a single point only

        \note Only for dim==0
     */
    AxisAlignedCubeGeometry(const Dune::FieldVector<ctype,coorddim> lower)
      : lower_(lower)
    {}

    /** \brief Type of the cube.  Here: a hypercube of the correct dimension */
    GeometryType type() const
    {
      return GeometryTypes::cube(dim);
    }

    /** \brief Map a point in local (element) coordinates to world coordinates */
    GlobalCoordinate global(const LocalCoordinate& local) const
    {
      GlobalCoordinate result;
      if (dim == coorddim) {        // fast case
        for (size_t i=0; i<coorddim; i++)
          result[i] = lower_[i] + local[i]*(upper_[i] - lower_[i]);
      } else if (dim == 0) {              // a vertex -- the other fast case
        result = lower_;          // hope for named-return-type-optimization
      } else {          // slow case
        size_t lc=0;
        for (size_t i=0; i<coorddim; i++)
          result[i] = (axes_[i])
                      ? lower_[i] + local[lc++]*(upper_[i] - lower_[i])
                      : lower_[i];
      }
      return result;
    }

    /** \brief Map a point in global (world) coordinates to element coordinates */
    LocalCoordinate local(const GlobalCoordinate& global) const
    {
      LocalCoordinate result;
      if (dim == coorddim) {        // fast case
        for (size_t i=0; i<dim; i++)
          result[i] = (global[i] - lower_[i]) / (upper_[i] - lower_[i]);
      } else if (dim != 0) {          // slow case
        size_t lc=0;
        for (size_t i=0; i<coorddim; i++)
          if (axes_[i])
            result[lc++] = (global[i] - lower_[i]) / (upper_[i] - lower_[i]);
      }
      return result;
    }

    /** \brief Jacobian transposed of the transformation from local to global coordinates */
    JacobianTransposed jacobianTransposed([[maybe_unused]] const LocalCoordinate& local) const
    {
      JacobianTransposed result;

      // Actually compute the result.  Uses different methods depending
      // on what kind of matrix JacobianTransposed is.
      jacobianTransposed(result);

      return result;
    }

    /** \brief Inverse Jacobian transposed of the transformation from local to global coordinates */
    JacobianInverseTransposed jacobianInverseTransposed([[maybe_unused]] const LocalCoordinate& local) const
    {
      JacobianInverseTransposed result;

      // Actually compute the result.  Uses different methods depending
      // on what kind of matrix JacobianTransposed is.
      jacobianInverseTransposed(result);

      return result;
    }

    /** \brief Jacobian of the transformation from local to global coordinates */
    Jacobian jacobian([[maybe_unused]] const LocalCoordinate& local) const
    {
      return jacobianTransposed(local).transposed();
    }

    /** \brief Inverse Jacobian of the transformation from local to global coordinates */
    JacobianInverse jacobianInverse([[maybe_unused]] const LocalCoordinate& local) const
    {
      return jacobianInverseTransposed(local).transposed();
    }

    /** \brief Return the integration element, i.e., the determinant term in the integral
               transformation formula
     */
    Volume integrationElement([[maybe_unused]] const LocalCoordinate& local) const
    {
      return volume();
    }

    /** \brief Return center of mass of the element */
    GlobalCoordinate center() const
    {
      GlobalCoordinate result;
      if (dim==0)
        result = lower_;
      else {
        // Since lower_==upper_ for unused coordinates, this always does the right thing
        for (size_t i=0; i<coorddim; i++)
          result[i] = CoordType(0.5) * (lower_[i] + upper_[i]);
      }
      return result;
    }

    /** \brief Return the number of corners of the element */
    int corners() const
    {
      return 1<<dim;
    }

    /** \brief Return world coordinates of the k-th corner of the element */
    GlobalCoordinate corner(int k) const
    {
      GlobalCoordinate result;
      if (dim == coorddim) {         // fast case
        for (size_t i=0; i<coorddim; i++)
          result[i] = (k & (1<<i)) ? upper_[i] : lower_[i];
      } else if (dim == 0) {         // vertex
        result = lower_;            // rely on named return-type optimization
      } else {                // slow case
        unsigned int mask = 1;

        for (size_t i=0; i<coorddim; i++) {
          if (not axes_[i])
            result[i] = lower_[i];
          else {
            result[i] = (k & mask) ? upper_[i] : lower_[i];
            mask = (mask<<1);
          }
        }
      }


      return result;
    }

    /** \brief Return the element volume */
    Volume volume() const
    {
      ctype vol = 1;
      if (dim == coorddim) {         // fast case
        for (size_t i=0; i<dim; i++)
          vol *= upper_[i] - lower_[i];
        // do nothing if dim == 0
      } else if (dim != 0) {         // slow case
        for (size_t i=0; i<coorddim; i++)
          if (axes_[i])
            vol *= upper_[i] - lower_[i];
      }
      return vol;
    }

    /** \brief Return if the element is affine.  Here: yes */
    bool affine() const
    {
      return true;
    }

    friend Dune::Transitional::ReferenceElement< ctype, Dim<dim> > referenceElement ( const AxisAlignedCubeGeometry & /* geometry */ )
    {
      return ReferenceElements< ctype, dim >::cube();
    }

  private:
    // jacobianTransposed: fast case --> diagonal matrix
    void jacobianTransposed ( DiagonalMatrix<ctype,dim> &jacobianTransposed ) const
    {
      for (size_t i=0; i<dim; i++)
        jacobianTransposed.diagonal()[i] = upper_[i] - lower_[i];
    }

    // jacobianTransposed: slow case --> dense matrix
    void jacobianTransposed ( FieldMatrix<ctype,dim,coorddim> &jacobianTransposed ) const
    {
      if (dim==0)
          return;

      size_t lc = 0;
      for (size_t i=0; i<coorddim; i++)
        if (axes_[i])
          jacobianTransposed[lc++][i] = upper_[i] - lower_[i];
    }

    // jacobianInverseTransposed: fast case --> diagonal matrix
    void jacobianInverseTransposed ( DiagonalMatrix<ctype,dim> &jacobianInverseTransposed ) const
    {
      for (size_t i=0; i<dim; i++)
        jacobianInverseTransposed.diagonal()[i] = CoordType(1.0) / (upper_[i] - lower_[i]);
    }

    // jacobianInverseTransposed: slow case --> dense matrix
    void jacobianInverseTransposed ( FieldMatrix<ctype,coorddim,dim> &jacobianInverseTransposed ) const
    {
      if (dim==0)
          return;

      size_t lc = 0;
      for (size_t i=0; i<coorddim; i++)
        if (axes_[i])
          jacobianInverseTransposed[i][lc++] = CoordType(1.0) / (upper_[i] - lower_[i]);
    }

    Dune::FieldVector<ctype,coorddim> lower_;

    Dune::FieldVector<ctype,coorddim> upper_;

    std::bitset<coorddim> axes_;
  };

} // namespace Dune
#endif
