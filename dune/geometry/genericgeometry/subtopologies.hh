// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_SUBTOPOLOGIES_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_SUBTOPOLOGIES_HH

#include <cassert>
#include <vector>

#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>
#include <dune/common/unused.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, unsigned int codim >
    struct Size;

    template< class Topology, unsigned int codim, unsigned int i >
    struct SubTopology;



    // Size
    // ----

    template< class Topology, unsigned int dim, unsigned int codim >
    class SizeImpl;

    template< unsigned int dim, unsigned int codim >
    class SizeImpl< Point, dim, codim >
    {
      typedef Point Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((codim <= dim), "Invalid codimension");

    public:
      enum { value = 1 };
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim >
    class SizeImpl< Prism< BaseTopology >, dim, codim >
    {
      typedef Prism< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((codim <= dim), "Invalid codimension");

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

    public:
      enum { value = n + 2*m };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Prism< BaseTopology >, dim, 0 >
    {
      typedef Prism< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");

    public:
      enum { value = 1 };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Prism< BaseTopology >, dim, dim >
    {
      typedef Prism< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");

      enum { m = Size< BaseTopology, dim-1 > :: value };

    public:
      enum { value = 2*m };
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim >
    class SizeImpl< Pyramid< BaseTopology >, dim, codim >
    {
      typedef Pyramid< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((codim <= dim), "Invalid codimension");

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

    public:
      enum { value = m+n };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Pyramid< BaseTopology >, dim, 0 >
    {
      typedef Pyramid< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");

    public:
      enum { value = 1 };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Pyramid< BaseTopology >, dim, dim >
    {
      typedef Pyramid< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");

      enum { m = Size< BaseTopology, dim-1 > :: value };

    public:
      enum { value = m+1 };
    };

    /** \brief Statically compute the number of subentities of a given codimension */
    template< class Topology, unsigned int codim >
    struct Size
    {
      enum { value = SizeImpl< Topology, Topology :: dimension, codim > :: value };
    };




    /** \brief Compute the number of subentities of a given codimension */
    unsigned int size ( unsigned int topologyId, int dim, int codim );



    // SubTopology
    // -----------

    template< class Topology, unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl;

    template< unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl< Point, dim, codim, i >
    {
      typedef Point Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((codim <= dim), "Invalid codimension");
      static_assert((i < Size< Topology, codim > :: value), "Invalid subentity index");

    public:
      typedef Topology type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl< Prism< BaseTopology >, dim, codim, i >
    {
      typedef Prism< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((codim <= dim), "Invalid codimension");
      static_assert((i < Size< Topology, codim > :: value), "Invalid subentity index");

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef Prism< typename SubTopology< BaseTopology, codim, i > :: type > type;
      };

      template< bool >
      struct BaseSub
      {
        typedef typename SubTopology< BaseTopology, codim-1, i-(n+s*m) > :: type type;
      };

    public:
      typedef typename std::conditional< (i < n), PrismSub<true>, BaseSub<false> > :: type :: type type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Prism< BaseTopology >, dim, 0, i >
    {
      typedef Prism< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((i < Size< Topology, 0 > :: value), "Invalid subentity index");
    public:
      typedef Topology type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Prism< BaseTopology >, dim, dim, i >
    {
      typedef Prism< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((i < Size< Topology, dim > :: value), "Invalid subentity index");
    public:
      typedef Point type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl< Pyramid< BaseTopology >, dim, codim, i >
    {
      typedef Pyramid< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((codim <= dim), "Invalid codimension" );
      static_assert((i < Size< Topology, codim > :: value), "Invalid subentity index");

      enum { m = Size< BaseTopology, codim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        typedef typename SubTopology< BaseTopology, codim-1, i > :: type type;
      };

      template< bool >
      struct PyramidSub
      {
        typedef Pyramid< typename SubTopology< BaseTopology, codim, i-m > :: type > type;
      };

    public:
      typedef typename std::conditional< (i < m), BaseSub<true>, PyramidSub<false> > :: type :: type type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Pyramid< BaseTopology >, dim, 0, i >
    {
      typedef Pyramid< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((i < Size< Topology, 0 > :: value), "Invalid subentity index");

    public:
      typedef Topology type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Pyramid< BaseTopology >, dim, dim, i >
    {
      typedef Pyramid< BaseTopology > Topology;
      static_assert((dim == Topology :: dimension), "Wrong dimension");
      static_assert((i < Size< Topology, dim > :: value), "Invalid subentity index");

    public:
      typedef Point type;
    };

    template< class Topology, unsigned int codim, unsigned int i >
    struct SubTopology
    {
      typedef typename SubTopologyImpl< Topology, Topology :: dimension, codim, i > :: type type;
    };



    /** \brief Compute the topology id of a given subentity
     *
     * \param topologyId Topology id of entity
     * \param dim Dimension of entity
     * \param codim Codimension of the subentity that we are interested in
     * \param i Number of the subentity that we are interested in
     */
    unsigned int subTopologyId ( unsigned int topologyId, int dim, int codim, unsigned int i );




    // subTopologyNumbering
    // --------------------

    void subTopologyNumbering ( unsigned int topologyId, int dim, int codim, unsigned int i, int subcodim,
                                unsigned int *beginOut, unsigned int *endOut );

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_SUBTOPOLOGIES_HH
