// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_STATICGEOMETRYTYPE_HH
#define DUNE_GEOMETRY_STATICGEOMETRYTYPE_HH

#include <dune/common/std/type_traits.hh>
#include <dune/common/indices.hh>
#include <dune/common/keywords.hh>
#include <dune/geometry/type.hh>



namespace Dune {

/**
 * \brief A geometry type template encoding the geometry type statically
 *
 * This works like GeometryType but encodes topology id, dimension,
 * and the none-flag statically as template parameters. For interoperability
 * on object of any instance of StaticGeometryType can be converted to
 * the corresponding dynamic GeometryType object.
 */
template<unsigned int topologyId_T, std::size_t dim_T, bool none_T>
class StaticGeometryType
{
  using NormalizedType = StaticGeometryType<((topologyId_T >> 1) << 1), dim_T, none_T>;
public:

  /** \brief Return the topology id of the type */
  static constexpr auto id()
  {
    return Dune::index_constant<topologyId_T>();
  }

  /** \brief Return dimension of the type */
  static constexpr auto dim()
  {
    return Dune::index_constant<dim_T>();
  }

  /** \brief Return true if entity is a singular of any dimension */
  static constexpr auto isNone()
  {
    return Dune::Std::bool_constant<none_T>();
  }


  /** \brief Return true if entity is a vertex */
  static constexpr auto isVertex()
  {
    return Dune::Std::bool_constant<dim()==0>();
  }

  /** \brief Return true if entity is a line segment */
  static constexpr auto isLine()
  {
    return Dune::Std::bool_constant<dim()==1>();
  }

  /** \brief Return true if entity is a triangle */
  static constexpr auto isTriangle()
  {
    return std::is_same<NormalizedType, StaticGeometryType<0, 2, false>>();
  }

  /** \brief Return true if entity is a quadrilateral */
  static constexpr auto isQuadrilateral()
  {
    return std::is_same<NormalizedType, StaticGeometryType<0b0010, 2, false>>();
  }

  /** \brief Return true if entity is a tetrahedron */
  static constexpr auto isTetrahedron()
  {
    return std::is_same<NormalizedType, StaticGeometryType<0, 3, false>>();
  }

  /** \brief Return true if entity is a pyramid */
  static constexpr auto isPyramid()
  {
    return std::is_same<NormalizedType, StaticGeometryType<0b0010, 3, false>>();
  }

  /** \brief Return true if entity is a prism */
  static constexpr auto isPrism()
  {
    return std::is_same<NormalizedType, StaticGeometryType<0b0100, 3, false>>();
  }

  /** \brief Return true if entity is a hexahedron */
  static constexpr auto isHexahedron()
  {
    return std::is_same<NormalizedType, StaticGeometryType<0b0110, 3, false>>();
  }

  /** \brief Return true if entity is a simplex of any dimension */
  static constexpr auto isSimplex()
  {
    return Dune::Std::bool_constant<(!isNone()) && ((id() | 1) == 1)>();
  }

  /** \brief Return true if entity is a cube of any dimension */
  static constexpr auto isCube()
  {
    return Dune::Std::bool_constant<(!isNone()) && ((id() ^ ((1 << dim())-1)) >> 1 == 0)>();
  }

  static constexpr Dune::GeometryType dynamicType()
  {
    return Dune::GeometryType(id(), dim(), isNone());
  }

  constexpr operator Dune::GeometryType () const
  {
    return dynamicType();
  }

};



namespace StaticGeometryTypes {


  // Simplex *******************************************************************
  template<std::size_t dim>
  using Simplex = StaticGeometryType<0,dim,false>;

  template<std::size_t dim>
  inline constexpr auto simplex()
  {
    return Simplex<dim>();
  }

  template<std::size_t dim>
  inline constexpr auto simplex(Dune::index_constant<dim>)
  {
    return Simplex<dim>();
  }

  // Cube **********************************************************************
  template<std::size_t dim>
  using Cube = StaticGeometryType<((dim>1) ? ((1 << dim) - 2) : 0),dim,false>;

  template<std::size_t dim>
  inline constexpr auto cube()
  {
    return Cube<dim>();
  }

  template<std::size_t dim>
  inline constexpr auto cube(Dune::index_constant<dim>)
  {
    return Cube<dim>();
  }

  // None **********************************************************************
  template<std::size_t dim>
  using None = StaticGeometryType<0, dim, true>;

  template<std::size_t dim>
  inline constexpr auto none()
  {
    return None<dim>();
  }

  template<std::size_t dim>
  inline constexpr auto none(Dune::index_constant<dim>)
  {
    return None<dim>();
  }

#ifndef __cpp_inline_variables
    namespace {
#endif

      using Vertex = StaticGeometryType<0,0,false>;
      DUNE_INLINE_VARIABLE constexpr auto vertex = Vertex();

      using Line = StaticGeometryType<0,1,false>;
      DUNE_INLINE_VARIABLE constexpr auto line = Line();

      using Triangle = StaticGeometryTypes::Simplex<2>;
      DUNE_INLINE_VARIABLE constexpr auto triangle = Triangle();

      using Quadrilateral = StaticGeometryTypes::Cube<2>;
      DUNE_INLINE_VARIABLE constexpr auto quadrilateral = Quadrilateral();

      using Tetrahedron = StaticGeometryTypes::Simplex<3>;
      DUNE_INLINE_VARIABLE constexpr auto tetrahedron = Tetrahedron();

      using Pyramid = StaticGeometryType<0b0010,3,false>;
      DUNE_INLINE_VARIABLE constexpr auto pyramid = Pyramid();

      using Prism = StaticGeometryType<0b0100,3,false>;
      DUNE_INLINE_VARIABLE constexpr auto prism = Prism();

      using Hexahedron = StaticGeometryTypes::Cube<3>;
      DUNE_INLINE_VARIABLE constexpr auto hexahedron = Hexahedron();

#ifndef __cpp_inline_variables
    }
#endif

} // namespace StaticGeometryTypes




//! Compute per-dimension indices for geometry types
class LocalHybridGeometryTypeIndex
{
  inline static constexpr std::size_t regular_size(std::size_t dim)
  {
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
    return regular_size(dim) + 1;
  }

  /**
   * \brief Compute the index for the given GeometryType
   */
  inline static constexpr std::size_t index(const GeometryType &gt)
  {
    return gt.isNone() ?  regular_size(gt.dim()) : (gt.id() >> 1);
  }

  /**
   * \brief Compute the index for the given GeometryType
   */
  template<unsigned int topologyId, std::size_t dim, bool none>
  inline static constexpr auto index(const StaticGeometryType<topologyId, dim, none>& gt)
  {
    return Dune::index_constant<none ?  regular_size(dim) : (topologyId >> 1)>();
  }

  //! compute the geometry type for the given local index and dimension
  inline static constexpr GeometryType type(std::size_t dim, std::size_t index) {
    return (index == regular_size(dim)) ?
      GeometryTypes::none(dim) :
      GeometryType(static_cast< unsigned int >(index << 1), dim);
  }

  //! compute the geometry type for the given local index and dimension
  template<std::size_t dim, std::size_t index>
  inline static constexpr auto type()
  {
    return std::conditional_t<(index == regular_size(dim)),
           StaticGeometryTypes::None<dim>,
           StaticGeometryType<(index << 1), dim, false>>();
  }

  template<std::size_t dim, std::size_t index>
  inline static constexpr auto type(Dune::index_constant<dim>, Dune::index_constant<index>)
  {
    return type<dim, index>();
  }

};



} // namespace Dune



#endif // DUNE_GEOMETRY_STATICGEOMETRYTYPE_HH
