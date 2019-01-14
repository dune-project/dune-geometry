// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GEOMETRY_IDENTITYGEOMETRY_HH
#define DUNE_GEOMETRY_IDENTITYGEOMETRY_HH

/** \file
    \brief A geometry implementation for identity mappings
 */

#include <dune/common/fvector.hh>
#include <dune/common/diagonalmatrix.hh>
#include <dune/common/unused.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>


namespace Dune {

  /** \brief A geometry implementation for identity mappings
   *
   * \tparam dim Dimension of the cube
   * \tparam coorddim Dimension of the space that the cube lives in
   */
  template <class CoordType, unsigned int dim>
  class IdentityGeometry
  {


  public:

    /** \brief Dimension of the cube element */
    static constexpr int mydimension = dim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr int coorddimension = dim;

    /** \brief Type used for single coordinate coefficients */
    using ctype = CoordType;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = FieldVector<ctype,dim>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = FieldVector<ctype,dim>;

    /** \brief Return type of jacobianTransposed

        This is a fast DiagonalMatrix if dim==coorddim, and a FieldMatrix otherwise.
        The FieldMatrix will never contain more than one entry per row,
        hence it could be replaced by something more efficient.
     */
    using JacobianTransposed = DiagonalMatrix<ctype,dim>;

    /** \brief Return type of jacobianInverseTransposed

        This is a fast DiagonalMatrix if dim==coorddim, and a FieldMatrix otherwise.
        The FieldMatrix will never contain more than one entry per column,
        hence it could be replaced by something more efficient.
     */
    using JacobianInverseTransposed = DiagonalMatrix<ctype,dim>;

    IdentityGeometry(GeometryType type)
      : _type(type)
    {}

    /** \brief Type of the cube.  Here: a hypercube of the correct dimension */
    GeometryType type() const
    {
      return _type;
    }

    /** \brief Map a point in local (element) coordinates to world coordinates */
    GlobalCoordinate global(const LocalCoordinate& local) const
    {
      return local;
    }

    /** \brief Map a point in global (world) coordinates to element coordinates */
    LocalCoordinate local(const GlobalCoordinate& global) const
    {
      return global;
    }

    /** \brief Jacobian transposed of the transformation from local to global coordinates */
    JacobianTransposed jacobianTransposed(DUNE_UNUSED const LocalCoordinate& local) const
    {
      return {ctype(1.0)};
    }

    /** \brief Jacobian transposed of the transformation from local to global coordinates */
    JacobianInverseTransposed jacobianInverseTransposed(DUNE_UNUSED const LocalCoordinate& local) const
    {
      return {ctype(1.0)};
    }

    /** \brief Return the integration element, i.e., the determinant term in the integral
               transformation formula
     */
    ctype integrationElement(DUNE_UNUSED const LocalCoordinate& local) const
    {
      return 1.0;
    }

    /** \brief Return center of mass of the element */
    GlobalCoordinate center() const
    {
      return referenceElement<ctype>(_type,Dim<dim>{}).position(0,0);
    }

    /** \brief Return the number of corners of the element */
    int corners() const
    {
      return referenceElement<ctype>(_type,Dim<dim>{}).size(dim);
    }

    /** \brief Return world coordinates of the k-th corner of the element */
    GlobalCoordinate corner(int k) const
    {
      return referenceElement<ctype>(_type,Dim<dim>{}).position(dim,k);
    }

    /** \brief Return the element volume */
    ctype volume() const
    {
      return referenceElement<ctype>(_type,Dim<dim>{}).volume();
    }

    /** \brief Return if the element is affine.  Here: yes */
    bool affine() const
    {
      return true;
    }

    friend Dune::Transitional::ReferenceElement< ctype, Dim<dim> > referenceElement ( const IdentityGeometry &geometry )
    {
      return referenceElement<ctype>(geometry._type,Dim<dim>{});
    }

  private:

    GeometryType _type;

  };

} // namespace Dune
#endif
