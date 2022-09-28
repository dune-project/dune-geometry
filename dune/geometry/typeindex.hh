// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TYPEINDEX_HH
#define DUNE_GEOMETRY_TYPEINDEX_HH

/**
 * \file
 * \brief Helper classes to provide indices for geometrytypes for use in a
 *        vector
 */

#include <cstddef>

#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include "type.hh"

namespace Dune
{
  //! Compute per-dimension indices for geometry types
  class LocalGeometryTypeIndex
  {
    /**
     * \brief Compute the number of regular geometry types for the given dimension
     *
     * Regular geometry type are those which have a topologyId, i.e. "None" is
     * not a regular geometry type.
     */
    inline static constexpr std::size_t regular_size(std::size_t dim)
    {
      // The following expression is derived from the expression for
      // GlobalGeometryTypeIndex::regular_offset().  Subtracting
      // regular_offset(dim+1)-regular_offset(dim) we get:
      //
      //   ((1 << dim+1) >> 1) - ((1 << dim) >> 1)
      //
      // We always have
      //
      //   dim >= 0,
      //
      // so
      //
      //   (1 << dim+1) >= 2   and   (1 << dim+2) % 2 == 0.
      //
      // So if we apply a single right-shift to that, we will never lose any
      // set bits, thus
      //
      //   ((1 << dim+1) >> 1) == (1 << dim)
      return (1 << dim) - ((1 << dim) >> 1);
    }

  public:
    /**
     * \brief Compute total number of geometry types for the given dimension
     *
     * This includes irregular geometry types such as "None".
     */
    inline static constexpr std::size_t size(std::size_t dim)
    {
      // one for "none"
      return regular_size(dim) + 1;
    }

    /**
     * \brief Compute the index for the given geometry type within its dimension
     *
     * Geometry types from different dimensions my get the same index. If that
     * is not what you want, maybe you should look at GlobalGeometryTypeIndex.
     */
    inline static constexpr std::size_t index(const GeometryType &gt)
    {
      return gt.isNone() ?  regular_size(gt.dim()) : (gt.id() >> 1);
    }

    //! compute the geometry type for the given local index and dimension
    inline static constexpr GeometryType type(std::size_t dim, std::size_t index) {
      return (index == regular_size(dim)) ?
        GeometryTypes::none(dim) :
        // the cast to unsigned makes sure this is interpreted as the topology
        // ID constructor
        GeometryType(static_cast< unsigned int >(index << 1), dim);
    }
  };

  //! Compute indices for geometry types, taking the dimension into account
  class GlobalGeometryTypeIndex
  {
    /**
     * \brief Compute the starting index for a given dimension ignoring
     *        irregular geometry types
     *
     * This ignores irregular geometry types so it is not useful in itself.
     * Have a look at offset() which does include the irregular geometry types.
     */
    inline static constexpr std::size_t regular_offset(std::size_t dim)
    {
      // The number of regular geometry types in a given dimension is
      // 2^(dim-1).  For dim==0 this would yield 1/2 geometry types (which is
      // obviously bogus, dim==0 has one regular geometry type, the point).
      // The following expression relies on 1 >> 1 == 0 to treat dim==0
      // specially.
      return (1 << dim) >> 1;
    }

  public:
    /**
     * \brief Compute the starting index for a given dimension including
     *        irregular geometry types
     */
    inline static constexpr std::size_t offset(std::size_t dim)
    {
      // dim times "none"
      return regular_offset(dim) + dim;
    }

    /**
     * \brief Compute total number of geometry types up to and including the
     *        given dimension
     *
     * This includes irregular geometry types such as "None".
     */
    inline static constexpr std::size_t size(std::size_t maxdim)
    {
      return offset(maxdim+1);
    }

    /**
     * \brief Compute the index for the given geometry type over all dimensions
     *
     * Geometry types from different dimensions will get different indices,
     * and lower dimensions will always have lower indices than higher
     * dimensions.  If that is not what you want, maybe you should look at
     * LocalGeometryTypeIndex.
     */
    inline static constexpr std::size_t index(const GeometryType &gt)
    {
      return offset(gt.dim()) + LocalGeometryTypeIndex::index(gt);
    }
  };

  namespace Impl {

    // Map a dynamic GeometryType to a static integral_constant<GeometryType::Id, ...>
    template<int dim, class F>
    static auto toGeometryTypeIdConstant(const GeometryType& gt, F&& f) {
      // Transform LocalGeometryTypeIndex to GeometryType::Id
      auto callWithId = [&](auto index) {
        static constexpr auto id = LocalGeometryTypeIndex::type(dim, decltype(index)::value).toId();
        return f(std::integral_constant<GeometryType::Id, id>{});
      };
      // switchCases needs a fallback to determine the return type.
      auto fallBack = [&]() { return callWithId(Dune::Indices::_0); };
      // Iterate over all _regular_ GeometryType indices for given dimension
      auto allIndices = std::make_index_sequence<LocalGeometryTypeIndex::size(dim)-1>{};
      return Dune::Hybrid::switchCases(allIndices, LocalGeometryTypeIndex::index(gt), callWithId, fallBack);
    }

  } // namespace Impl
} // namespace Dune

#endif // DUNE_GEOMETRY_TYPEINDEX_HH
