// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTS_HH

#include <cassert>

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/array.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelementimplementation.hh>

namespace Dune
{

  // ReferenceElementContainer
  // -------------------------

  template< class ctype, int dim >
  class ReferenceElementContainer
  {
    static const unsigned int numTopologies = (1u << dim);

  public:
    typedef ReferenceElement< ctype, dim > value_type;
    typedef const value_type *const_iterator;

    ReferenceElementContainer ()
    {
      for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        values_[ topologyId ].initialize( topologyId );
    }

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == dim );
      return values_[ type.id() ];
    }

    const value_type &simplex () const
    {
      return values_[ Impl::SimplexTopology< dim >::type::id ];
    }

    const value_type &cube () const
    {
      return values_[ Impl::CubeTopology< dim >::type::id ];
    }

    const value_type &pyramid () const
    {
      return values_[ Impl::PyramidTopology< dim >::type::id ];
    }

    const value_type &prism () const
    {
      return values_[ Impl::PrismTopology< dim >::type::id ];
    }

    const_iterator begin () const { return values_; }
    const_iterator end () const { return values_ + numTopologies; }

  private:
    value_type values_[ numTopologies ];
  };



  // ReferenceElements
  // ------------------------

  /** \brief Class providing access to the singletons of the
   *  reference elements
   *
   *  Special methods are available for
   *  simplex and cube elements of any dimension.
   *  The method general can be used to obtain the reference element
   *  for a given geometry type.
   *
   *  \ingroup GeometryReferenceElements
   */
  template< class ctype, int dim >
  struct ReferenceElements
  {
    typedef typename ReferenceElementContainer< ctype, dim >::const_iterator Iterator;

    //! get general reference elements
    static const ReferenceElement< ctype, dim > &
    general ( const GeometryType &type )
    {
      return container() ( type );
    }

    //! get simplex reference elements
    static const ReferenceElement< ctype, dim > &simplex ()
    {
      return container().simplex();
    }

    //! get hypercube reference elements
    static const ReferenceElement< ctype, dim > &cube ()
    {
      return container().cube();
    }

    static Iterator begin () { return container().begin(); }
    static Iterator end () { return container().end(); }

  private:
    DUNE_EXPORT static const ReferenceElementContainer< ctype, dim > &container ()
    {
      static ReferenceElementContainer< ctype, dim > container;
      return container;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
