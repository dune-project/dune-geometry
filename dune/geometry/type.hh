// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_TYPE_HH
#define DUNE_GEOMETRY_TYPE_HH

/** \file
 *  \brief A unique label for each type of element that can occur in a grid
 */

#include <cassert>

#include <string>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/keywords.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/unused.hh>

namespace Dune
{

  // forward declaration needed for deprecated makeFromVertices
  class GeometryType;
  GeometryType geometryTypeFromVertexCount(unsigned int dim, unsigned int vertices);

  namespace Impl
  {

    enum TopologyConstruction { pyramidConstruction = 0, prismConstruction = 1 };



    // Basic Topology Types
    // --------------------

    struct Point
    {
      static const unsigned int dimension = 0;
      static const unsigned int numCorners = 1;

      static const unsigned int id = 0;

      static std::string name () { return "p"; }
    };


    template< class BaseTopology >
    struct Prism
    {
      static const unsigned int dimension = BaseTopology::dimension + 1;
      static const unsigned int numCorners = 2 * BaseTopology::numCorners;

      static const unsigned int id = BaseTopology::id | ((unsigned int)prismConstruction << (dimension-1));

      static std::string name () { return BaseTopology::name() + "l"; }
    };


    template< class BaseTopology >
    struct Pyramid
    {
      static const unsigned int dimension = BaseTopology::dimension + 1;
      static const unsigned int numCorners = BaseTopology::numCorners + 1;

      static const unsigned int id = BaseTopology::id | ((unsigned int)pyramidConstruction << (dimension-1));

      static std::string name () { return BaseTopology::name() + "o"; }
    };



    // Properties of Topologies
    // ------------------------

    template< class Topology >
    struct IsSimplex
      : public std::integral_constant< bool, (Topology::id >> 1) == 0 >
    {};

    template< class Topology >
    struct IsCube
      : public std::integral_constant< bool,  (Topology::id | 1) == (1 << Topology::dimension) - 1 >
    {};



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

    /** \brief check whether a specific topology construction was used to create a
     *         given codimension
     *
     *  \param[in]  construction  construction to check for
     *  \param[in]  topologyId    id of the topology
     *  \param[in]  dim           dimension of the topology
     *  \param[in]  codim         codimension for which the information is desired
     *                            (defaults to 0)
     *
     *  \returns true, if construction was used to generate the codimension the
     *           topology.
     */
    inline static bool isTopology ( TopologyConstruction construction, unsigned int topologyId, int dim, int codim = 0 ) noexcept
    {
      assert( (dim > 0) && (topologyId < numTopologies( dim )) );
      assert( (0 <= codim) && (codim <= dim) );
      return (codim >= (dim-1)) || (((topologyId >> (dim-codim-1)) & 1) == (unsigned int)construction);
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



    // SimplexTopology
    // ---------------

    template< unsigned int dim >
    struct SimplexTopology
    {
      typedef Pyramid< typename SimplexTopology< dim-1 >::type > type;
    };

    template<>
    struct SimplexTopology< 0 >
    {
      typedef Point type;
    };



    // CubeTopology
    // ------------

    template< unsigned int dim >
    struct CubeTopology
    {
      typedef Prism< typename CubeTopology< dim-1 >::type > type;
    };

    template<>
    struct CubeTopology< 0 >
    {
      typedef Point type;
    };



    // PyramidTopology
    // ---------------

    template< unsigned int dim >
    struct PyramidTopology
    {
      typedef Pyramid< typename CubeTopology< dim-1 >::type > type;
    };



    // PrismTopology
    // -------------

    template< unsigned int dim >
    struct PrismTopology
    {
      typedef Prism< typename SimplexTopology< dim-1 >::type > type;
    };




    // IfTopology
    // ----------

    template< template< class > class Operation, int dim, class Topology = Point >
    struct IfTopology
    {
      template< class... Args >
      static auto apply ( unsigned int topologyId, Args &&... args )
      {
        if( topologyId & 1 )
          return IfTopology< Operation, dim-1, Prism< Topology > >::apply( topologyId >> 1, std::forward< Args >( args )... );
        else
          return IfTopology< Operation, dim-1, Pyramid< Topology > >::apply( topologyId >> 1, std::forward< Args >( args )... );
      }
    };

    template< template< class > class Operation, class Topology >
    struct IfTopology< Operation, 0, Topology >
    {
      template< class... Args >
      static auto apply ( unsigned int topologyId, Args &&... args )
      {
        DUNE_UNUSED_PARAMETER( topologyId );
        return Operation< Topology >::apply( std::forward< Args >( args )... );
      }
    };

  } // namespace Impl



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

    /** \brief Topology Id element */
    unsigned int topologyId_;

    /** \brief Dimension of the element */
    unsigned char dim_  : 7;

    /** \brief bool if this is none-type */
    bool none_ : 1;

  public:

    /** @name Constructors */
    /*@{*/

    /** \brief Default constructor, not initializing anything */
    constexpr GeometryType ()
      : topologyId_(0), dim_(0), none_(true)
    {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    /** \brief Constructor, using the basic type and the dimension */
    GeometryType(BasicType basicType, unsigned int dim)
      DUNE_DEPRECATED_MSG("The GeometryType constructor taking BasicType is deprecated and will be removed after DUNE 2.6")
      : topologyId_(0), dim_(dim), none_((basicType == GeometryType::none) ? true : false)
    {
      if (dim < 2)
        return;
      switch( basicType )
      {
      case GeometryType::simplex :
        topologyId_ = 0;
        break;
      case GeometryType::cube :
        topologyId_ = ((1 << dim) - 1);
        break;
      case GeometryType::pyramid :
        if (dim == 3)
          topologyId_ = 0b0011;
        else
          DUNE_THROW( RangeError,
                      "Invalid basic geometry type: no pyramids for dimension " << dim << "." );
        break;
      case GeometryType::prism :
        if (dim == 3)
          topologyId_ = 0b0101;
        else
          DUNE_THROW( RangeError,
                      "Invalid basic geometry type: no prisms for dimension " << dim << "." );
        break;
      case GeometryType::none :
        break;
      default :
        DUNE_THROW( RangeError,
                    "Invalid basic geometry type: " << basicType << " for dimension " << dim << "." );
      }
    }
#pragma GCC diagnostic pop

    /** \brief Constructor, using the topologyId (integer), the dimension and a flag for type none.
     * \note With this constructor, you can easily create an invalid GeometryType,
     *       it is mostly here for internal use!
     *       the TypologyType, users are encouraged to use the
     *       GeometryType(TopologyType t) constructor.
     */
    constexpr GeometryType(unsigned int topologyId, unsigned int dim, bool none)
      : topologyId_(topologyId), dim_(dim), none_(none)
    {}

    /** \brief Constructor, using the topologyId (integer) and the dimension
     * \note the topologyId is a binary encoded representation of
     *       the TypologyType, users are encouraged to use the
     *       GeometryType(TopologyType t) constructor.
     */
    constexpr GeometryType(unsigned int topologyId, unsigned int dim)
      : topologyId_(topologyId), dim_(dim), none_(false)
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
      class = Dune::void_t<decltype(TopologyType::dimension), decltype(TopologyType::id)>>
    explicit GeometryType(TopologyType t)
      : topologyId_(TopologyType::id), dim_(TopologyType::dimension), none_(false)
    {
      DUNE_UNUSED_PARAMETER(t);
    }

    /** \brief Constructor for vertices and segments */
    explicit GeometryType(unsigned int dim)
      : topologyId_(0), dim_(dim), none_(false)
    {
      assert(dim < 2);
    }

    /** \brief Constructor for vertices and segments */
    // We need this constructor for "int" and "unsigned int",
    // because otherwise GeometryType(int) would try to call the
    // generic GeometryType(TopologyType) constructor
    explicit GeometryType(int dim)
      : topologyId_(0), dim_(dim), none_(false)
    {
      assert(dim < 2);
    }

    /** @} */


    /** @name Setup Methods */
    /*@{*/

    /** \brief Make a vertex */
    DUNE_DEPRECATED_MSG("makeVertex() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::vertex instead")
    void makeVertex() {
      none_  = false;
      dim_ = 0;
      topologyId_ = 0;
    }

    /** \brief Make a line segment */
    DUNE_DEPRECATED_MSG("makeLine() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::line instead")
    void makeLine() {
      none_  = false;
      dim_ = 1;
      topologyId_ = 0;
    }

    /** \brief Make a triangle */
    DUNE_DEPRECATED_MSG("makeTriangle() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::triangle instead")
    void makeTriangle() {
      none_  = false;
      dim_ = 2;
      topologyId_ = 0;
    }

    /** \brief Make a quadrilateral */
    DUNE_DEPRECATED_MSG("makeQuadrilateral() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::quadrilateral instead")
    void makeQuadrilateral() {
      none_  = false;
      dim_ = 2;
      topologyId_ = 0b0011;
    }

    /** \brief Make a tetrahedron */
    DUNE_DEPRECATED_MSG("makeTetrahedron() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::tetrahedron instead")
    void makeTetrahedron() {
      none_  = false;
      dim_ = 3;
      topologyId_ = 0;
    }

    /** \brief Make a pyramid */
    DUNE_DEPRECATED_MSG("makePyramid() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::pyramid instead")
    void makePyramid() {
      none_  = false;
      dim_ = 3;
      topologyId_ = 0b0011;
    }

    /** \brief Make a prism */
    DUNE_DEPRECATED_MSG("makePrism() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::prism instead")
    void makePrism() {
      none_  = false;
      dim_ = 3;
      topologyId_ = 0b0101;       // (1 << (dim_-1)) - 1;
    }

    /** \brief Make a hexahedron */
    DUNE_DEPRECATED_MSG("makeHexahedron() is deprecated in DUNE 2.6, please use Dune::GeometryTypes::hexahedron instead")
    void makeHexahedron() {
      none_  = false;
      dim_ = 3;
      topologyId_ = 0b0111;
    }

    /** \brief Make a simplex of given dimension */
    DUNE_DEPRECATED_MSG("makeSimplex(dim) is deprecated in DUNE 2.6, please use Dune::GeometryTypes::simplex(dim) instead")
    void makeSimplex(unsigned int dim) {
      none_  = false;
      dim_ = dim;
      topologyId_ = 0;
    }

    /** \brief Make a hypercube of given dimension */
    DUNE_DEPRECATED_MSG("makeCube(dim) is deprecated in DUNE 2.6, please use Dune::GeometryTypes::cube(dim) instead")
    void makeCube(unsigned int dim) {
      none_  = false;
      dim_ = dim;
      topologyId_ = ((dim>1) ? ((1 << dim) - 1) : 0);
    }

    /** \brief Make a singular of given dimension */
    DUNE_DEPRECATED_MSG("makeNone(dim) is deprecated in DUNE 2.6, please use Dune::GeometryTypes::none(dim) instead")
    void makeNone(unsigned int dim) {
      none_ = true;
      dim_ = dim;
      topologyId_  = 0;
    }

    /** \brief Construct the correct geometry type given the dimension and the number of vertices
     *  \note This code only works up to dimension 3.
     *        In higher dimensions the number of vertices does not uniquely identify the type of polyhedron.
     */
    void makeFromVertices(unsigned int dim, unsigned int vertices) DUNE_DEPRECATED_MSG("Use the utility function geometryTypeFromVertexCount(...) instead.")
    {
      *this = geometryTypeFromVertexCount(dim, vertices);
      return;
    }

    /*@}*/


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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  /** \brief Prints a GeometryType::BasicType to an output stream */
  inline std::ostream& operator<< (std::ostream& s, GeometryType::BasicType type)
  {
    switch (type) {
    case GeometryType::simplex :
      s << "simplex";
      break;
    case GeometryType::cube :
      s << "cube";
      break;
    case GeometryType::pyramid :
      s << "pyramid";
      break;
    case GeometryType::prism :
      s << "prism";
      break;
    case GeometryType::extended :
      s << "other";
    case GeometryType::none :
      s << "none";
      break;
    default :
      DUNE_THROW(Exception, "invalid GeometryType::BasicType");
    }
    return s;
  }
#pragma GCC diagnostic pop



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


} // namespace Dune

// include utility header needed for deprecated makeFromVertices
#include "utility/typefromvertexcount.hh"

#endif // DUNE_GEOMETRY_TYPE_HH
