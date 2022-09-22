// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TYPE_HH
#define DUNE_GEOMETRY_TYPE_HH

/** \file
 *  \brief A unique label for each type of element that can occur in a grid
 */

#include <cassert>

#include <string>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/keywords.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/unused.hh>

namespace Dune
{

  namespace Impl
  {

    enum TopologyConstruction { pyramidConstruction = 0, prismConstruction = 1 };

    // Dynamic Topology Properties
    // ---------------------------

    /** \brief obtain the number of topologies of a given dimension
     *
     *  \note Valid topology ids are 0,...,numTopologies(dim)-1.
     *
     *  \param[in]  dim  dimension
     *
     *  \returns number of topologies for the dimension
     */
    inline static unsigned int numTopologies ( int dim ) noexcept
    {
      return (1u << dim);
    }

    /** \brief check whether a pyramid construction was used to create a given
     *         codimension
     *
     *  \param[in]  topologyId  id of the topology
     *  \param[in]  dim         dimension of the topology
     *  \param[in]  codim       codimension for which the information is desired
     *                          (defaults to 0)
     *
     *  \returns true, if a pyramid construction was used to generate the
     *           codimension the topology.
     */
    inline bool static isPyramid ( unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim < dim) );
      return (((topologyId & ~1) & (1u << (dim-codim-1))) == 0);
    }

    /** \brief check whether a prism construction was used to create a given
     *         codimension
     *
     *  \param[in]  topologyId  id of the topology
     *  \param[in]  dim         dimension of the topology
     *  \param[in]  codim       codimension for which the information is desired
     *                          (defaults to 0)
     *
     *  \returns true, if a prism construction was used to generate the
     *           codimension the topology.
     */
    inline static bool isPrism ( unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim < dim) );
      return (( (topologyId | 1) & (1u << (dim-codim-1))) != 0);
    }

    /** \brief obtain the base topology of a given codimension
     *
     *  \param[in]  topologyId    id of the topology
     *  \param[in]  dim           dimension of the topology
     *  \param[in]  codim         codimension for which the information is desired
     *                            (defaults to 1)
     */
    inline static unsigned int baseTopologyId ( unsigned int topologyId, int dim, int codim = 1 ) noexcept
    {
      assert( (dim >= 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim <= dim) );
      return topologyId & ((1u << (dim-codim)) - 1);
    }

  } // namespace Impl

// the Topology classes are deprecated and will be removed for the 2.8.
// Temporarily a header 'deprecated_topology.hh' is provided which will be removed after the 2.9 release.
#if __GNUC__ >= 7
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#include <dune/geometry/deprecated_topology.hh>
#if __GNUC__ >= 7
#  pragma GCC diagnostic pop
#endif

  // GeometryType
  // -------------

  /** \brief Unique label for each type of entities that can occur in DUNE grids
   *
   * This class has to be extended if a grid implementation with new entity types
   * is added to DUNE.
   *
   * GeometryType is a C++ "literal type" and can be used in `constexpr` context if created
   * with a `constexpr` constructor.
   *
   * If you want to use a GeometryType as a template parameter, see GeometryType::Id.
   *
   * \ingroup GeometryType
   */
  class GeometryType
  {
  public:

    /** \brief Each entity can be tagged by one of these basic types
     *  plus its space dimension */
    enum
    BasicType {
      simplex,       //!< Simplicial element in any nonnegative dimension
      cube,          //!< Cube element in any nonnegative dimension
      pyramid,       //!< Four sided pyramid in three dimensions
      prism,         //!< Prism element in three dimensions
      extended,      //!< Other, more general topology, representable as topologyId
      none           //!< Even more general topology, cannot be specified by a topologyId. Two GeometryTypes with 'none' type are equal if and only if they have the same dimension.
    };

  private:

    /** \brief Dimension of the element */
    unsigned char dim_;

    /** \brief bool if this is none-type */
    bool none_;

    /** \brief Topology Id element */
    unsigned int topologyId_;

    // Internal type used for the Id. The exact nature of this type is kept
    // as an implementation detail on purpose. We use a scoped enum here because scoped enums
    // can be used as template parameters, but are not implicitly converted to other integral
    // types by the compiler. That way, we avoid unfortunate implicit conversion chains, e.g.
    // people trying to work with GlobalGeometryTypeIndex, but forgetting to actually call
    // GlobalGeometryTypeIndex::index(gt) and just using gt directly.
    enum class IdType : std::uint64_t
    {};

  public:

    /** \brief An integral id representing a GeometryType. */
    /**
     * Id is an unspecified built-in integral type that uniquely represents a GeometryType.
     * It mostly exists to be able to use a geometry type as a template parameter, as C++
     * does not let us use GeometryType directly for this purpose.
     *
     * GeometryType and GeometryType::Id are implicitly convertible to each other, while the
     * Id does not implicitly convert into other integral types. They should be used as follows:
     *
       \code
       // define a template with a GeometryType::Id parameter
       template<GeometryType::Id gtid>
       class Foo
       {
         // reconstruct a full-blown constexpr GeometryType as needed to access
         // information like the dimension etc.
         static constexpr GeometryType gt = gtid;
       };

       // Instantiate a Foo template
       Foo<GeometryTypes::triangle> foo;
       \endcode
     *
     * As you can see, the conversion between GeometryType and the id is completely transparent
     * to the user (apart from the slightly different template parameter type).
     *
     * \note The Id really only exists for this template parameter workaround. Do not use it to
     *       store a more compact version of the GeometryType - GeometryType and GeometryType::Id
     *       use the same amount of storage (64 bits).
     */
    using Id = IdType;

    /** \brief Construct an Id representing this GeometryType. */
    /**
     * This constructor exists mostly to transparently support using a GeometryType as a
     * template parameter.
     *
     * \sa Id
     */
    constexpr operator Id() const
    {
      // recreate the exact storage layout that this class is using, making conversion
      // extremely cheap
      std::uint64_t id = dim_ | (std::uint64_t(none_) << 8) | (std::uint64_t(topologyId_) << 32);
      return static_cast<Id>(id);
    }

    /** \brief Create an Id representation of this GeometryType. */
    /**
     * The returned Id encapsulates the whole information of this
     * GeometryType into an enum suitable for being used as template
     * parameter. The GeometryType can be reconstructed from the Id
     * using GeometryType{id}.
     *
     * This function was mainly introduced to support older GCC versions (<10.2).
     * There the implicit conversion from GeometryType to Id failed if a pure r-value
     * template argument based on a static class member was used.
     * (See dune/geometry/test/test-geometrytype-id.cc)
     */
    constexpr Id toId() const
    {
      return static_cast<Id>(*this);
    }

    /** \brief Reconstruct a Geometry type from a GeometryType::Id */
    /**
     * This constructor exists mostly to transparently support using a GeometryType as a
     * template parameter.
     *
     * \sa Id
     */
    constexpr GeometryType(Id id)
      : dim_(static_cast<std::uint64_t>(id) & 0xFF)
      , none_(static_cast<std::uint64_t>(id) & 0x100)
      , topologyId_(static_cast<std::uint64_t>(id) >> 32)
    {}

    /** @name Constructors */
    /*@{*/

    /** \brief Default constructor, not initializing anything */
    constexpr GeometryType ()
      : dim_(0), none_(true), topologyId_(0)
    {}

    /** \brief Constructor, using the topologyId (integer), the dimension and a flag for type none.
     * \note With this constructor, you can easily create an invalid GeometryType,
     *       it is mostly here for internal use!
     *       the TypologyType, users are encouraged to use the
     *       GeometryType(TopologyType t) constructor.
     */
    constexpr GeometryType(unsigned int topologyId, unsigned int dim, bool isNone)
      : dim_(dim), none_(isNone), topologyId_(topologyId)
    {}

    /** \brief Constructor, using the topologyId (integer) and the dimension
     * \note the topologyId is a binary encoded representation of
     *       the TypologyType, users are encouraged to use the
     *       GeometryType(TopologyType t) constructor.
     */
    constexpr GeometryType(unsigned int topologyId, unsigned int dim)
      : dim_(dim), none_(false), topologyId_(topologyId)
    {}

    /** \brief Constructor from static TopologyType class
     *
     * Constructs the GeometryType object from a static topology representation.
     *
     * \tparam TopologyType A class providing public static unsigned int members
     *                      TopologyType::dimension and TopologyType::id.
     *                      You can e.g. use the Point, Prism and Pyramid structs
     *                      from the Impl namespace.
     * \param t             Any object of type TopologyType. The object t itself is ignored.
     */
    template<class TopologyType,
      class = std::void_t<decltype(TopologyType::dimension), decltype(TopologyType::id)>>
    explicit GeometryType(TopologyType t)
      : dim_(TopologyType::dimension), none_(false), topologyId_(TopologyType::id)
    {
      DUNE_UNUSED_PARAMETER(t);
    }

    /** @} */


    /** @name Query Methods */
    /*@{*/
    /** \brief Return true if entity is a vertex */
    constexpr bool isVertex() const {
      return dim_==0;
    }

    /** \brief Return true if entity is a line segment */
    constexpr bool isLine() const {
      return dim_==1;
    }

    /** \brief Return true if entity is a triangle */
    constexpr bool isTriangle() const {
      return ! none_ && dim_==2 && (topologyId_ | 1) == 0b0001;
    }

    /** \brief Return true if entity is a quadrilateral */
    constexpr bool isQuadrilateral() const {
      return ! none_ && dim_==2 && (topologyId_ | 1) == 0b0011;
    }

    /** \brief Return true if entity is a tetrahedron */
    constexpr bool isTetrahedron() const {
      return ! none_ && dim_==3 && (topologyId_ | 1) == 0b0001;
    }

    /** \brief Return true if entity is a pyramid */
    constexpr bool isPyramid() const {
      return ! none_ && dim_==3 && (topologyId_ | 1) == 0b0011;
    }

    /** \brief Return true if entity is a prism */
    constexpr bool isPrism() const {
      return ! none_ && dim_==3 && (topologyId_ | 1) == 0b0101;
    }

    /** \brief Return true if entity is a hexahedron */
    constexpr bool isHexahedron() const {
      return ! none_ && dim_==3 && (topologyId_ | 1) == 0b0111;
    }

    /** \brief Return true if entity is a simplex of any dimension */
    constexpr bool isSimplex() const {
      return ! none_ && (topologyId_ | 1) == 1;
    }

    /** \brief Return true if entity is a cube of any dimension */
    constexpr bool isCube() const {
      return ! none_ && ((topologyId_ ^ ((1 << dim_)-1)) >> 1 == 0);
    }

    /** \brief Return true if entity was constructed with a conical product in the last step */
    constexpr bool isConical() const {
      return ! none_ && (((topologyId_ & ~1) & (1u << (dim_-1))) == 0);
    }

    /** \brief Return true if entity was constructed with a conical product in the chosen step
     *
     * \param step    0 <= step <= dim-1
     */
    constexpr bool isConical(const int& step) const {
      return ! none_ && (((topologyId_ & ~1) & (1u << step)) == 0);
    }

    /** \brief Return true if entity was constructed with a prismatic product in the last step */
    constexpr bool isPrismatic() const {
      return ! none_ && (( (topologyId_ | 1) & (1u << (dim_-1))) != 0);
    }

    /** \brief Return true if entity was constructed with a prismatic product in the chosen step
     *
     * \param step    0 <= step <= dim-1
     */
    constexpr bool isPrismatic(const int& step) const {
      return ! none_ && (( (topologyId_ | 1) & (1u << step)) != 0);
    }

    /** \brief Return true if entity is a singular of any dimension */
    constexpr bool isNone() const {
      return none_;
    }

    /** \brief Return dimension of the type */
    constexpr unsigned int dim() const {
      return dim_;
    }

    /** \brief Return the topology id of the type */
    constexpr unsigned int id() const {
      return topologyId_;
    }

    /*@}*/

    /** @name Comparison operators */

    /** \brief Check for equality. This method knows that in dimension 0 and 1
     *  all BasicTypes are equal.
     */
    constexpr bool operator==(const GeometryType& other) const {
      return ( ( none_ == other.none_ )
               && ( ( none_ == true )
                    || ( ( dim_ == other.dim_ )
                         && ( (topologyId_ >> 1) == (other.topologyId_ >> 1) )
                         )
                    )
               );
    }

    /** \brief Check for inequality */
    constexpr bool operator!=(const GeometryType& other) const {
      return ! ((*this)==other);
    }

    /** \brief less-than operation for use with maps */
    constexpr bool operator < (const GeometryType& other) const {
      return ( ( none_ < other.none_ )
               || ( !( other.none_ < none_ )
                    && ( ( dim_ < other.dim_ )
                         || ( (other.dim_ == dim_)
                              && ((topologyId_ >> 1) < (other.topologyId_ >> 1) )
                              )
                         )
                    )
               );
    }

    /*@}*/

  };

  /** \brief Prints the type to an output stream */
  inline std::ostream& operator<< (std::ostream& s, const GeometryType& a)
  {
    if (a.isSimplex())
    {
      s << "(simplex, " << a.dim() << ")";
      return s;
    }
    if (a.isCube())
    {
      s << "(cube, " << a.dim() << ")";
      return s;
    }
    if (a.isPyramid())
    {
      s << "(pyramid, 3)";
      return s;
    }
    if (a.isPrism())
    {
      s << "(prism, 3)";
      return s;
    }
    if (a.isNone())
    {
      s << "(none, " << a.dim() << ")";
      return s;
    }
    s << "(other [" << a.id() << "], " << a.dim() << ")";
    return s;
  }


  //! Predefined GeometryTypes for common geometries
  /**
   * \ingroup GeometryType
   * \related GeometryType
   */
  namespace GeometryTypes {

    //! Returns a GeometryType representing a simplex of dimension `dim`.
      /**
       * \ingroup GeometryType
       */
    inline constexpr GeometryType simplex(unsigned int dim)
    {
      return GeometryType(0,dim,false);
    }

    //! Returns a GeometryType representing a hypercube of dimension `dim`.
      /**
       * \ingroup GeometryType
       */
    inline constexpr GeometryType cube(unsigned int dim)
    {
      return GeometryType(((dim>1) ? ((1 << dim) - 1) : 0),dim,false);
    }

    //! Returns a GeometryType representing a singular of dimension `dim`.
      /**
       * \ingroup GeometryType
       */
    inline constexpr GeometryType none(unsigned int dim)
    {
      return GeometryType(0,dim,true);
    }

    /** \brief Return GeometryType of a conical construction with gt as base  */
    inline constexpr GeometryType conicalExtension(const GeometryType& gt)
    {
      return GeometryType(gt.id(), gt.dim()+1, gt.isNone());
    }

    /** \brief Return GeometryType of a prismatic construction with gt as base  */
    inline constexpr GeometryType prismaticExtension(const GeometryType& gt)
    {
      return GeometryType(gt.id() | ((1 << gt.dim())), gt.dim()+1, gt.isNone());
    }

#ifndef __cpp_inline_variables
    namespace {
#endif

      //! GeometryType representing a vertex.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType vertex = GeometryType(0,0,false);

      //! GeometryType representing a line.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType line = GeometryType(0,1,false);

      //! GeometryType representing a triangle.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType triangle = simplex(2);

      //! GeometryType representing a quadrilateral (a square).
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType quadrilateral = cube(2);

      //! GeometryType representing a tetrahedron.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType tetrahedron = simplex(3);

      //! GeometryType representing a 3D pyramid.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType pyramid = GeometryType(0b0011,3,false);

      //! GeometryType representing a 3D prism.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType prism = GeometryType(0b0101,3,false);

      //! GeometryType representing a hexahedron.
      /**
       * \ingroup GeometryType
       */
      DUNE_INLINE_VARIABLE constexpr GeometryType hexahedron = cube(3);

#ifndef __cpp_inline_variables
    }
#endif

  }

  namespace Impl
  {

    /** \brief Removes the bit for the highest dimension and returns the lower-dimensional GeometryType */
    inline constexpr GeometryType getBase(const GeometryType& gt) {
      return GeometryType(gt.id() & ((1 << (gt.dim()-1))-1), gt.dim()-1, gt.isNone());
    }


    // IfGeometryType
    // ----------

    template< template< GeometryType::Id > class Operation, int dim, GeometryType::Id geometryId = GeometryTypes::vertex >
    struct IfGeometryType
    {
      static constexpr GeometryType geometry = geometryId;
      template< class... Args >
      static auto apply ( GeometryType gt, Args &&... args )
      {
        GeometryType lowerGeometry(gt.id() >>1 , gt.dim()-1, gt.isNone());

        if( gt.id() & 1 )
          return IfGeometryType< Operation, dim-1, GeometryTypes::prismaticExtension(geometry).toId() >::apply( lowerGeometry, std::forward< Args >( args )... );
        else
          return IfGeometryType< Operation, dim-1, GeometryTypes::conicalExtension(geometry).toId() >::apply( lowerGeometry, std::forward< Args >( args )... );
      }
    };

    template< template< GeometryType::Id > class Operation, GeometryType::Id geometryId >
    struct IfGeometryType< Operation, 0, geometryId>
    {
      template< class... Args >
      static auto apply ([[maybe_unused]] GeometryType gt, Args &&... args )
      {
        return Operation< geometryId >::apply( std::forward< Args >( args )... );
      }
    };
  } // namespace Impl
} // namespace Dune

#endif // DUNE_GEOMETRY_TYPE_HH
