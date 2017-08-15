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
#include <array>

#include <dune/common/visibility.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/referenceelementimplementation.hh>

namespace Dune
{

  namespace Geo
  {

#ifndef DOXYGEN

    namespace Impl
    {

      // ReferenceElementContainer
      // -------------------------

      template< class ctype, int dim >
      class ReferenceElementContainer
      {
        static const unsigned int numTopologies = (1u << dim);

        using Implementation   = ReferenceElementImplementation< ctype, dim >;

      public:

        using ReferenceElement = Dune::Geo::ReferenceElement< Implementation >;
        using value_type       = ReferenceElement;
        using const_iterator   = const value_type*;

        ReferenceElementContainer ()
        {
          for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
            {
              implementations_[ topologyId ].initialize( topologyId );
              reference_elements_[ topologyId ].setImplementation( implementations_[ topologyId ] );
            }
        }

        const ReferenceElement& operator() ( const GeometryType &type ) const
        {
          assert( type.dim() == dim );
          return reference_elements_[ type.id() ];
        }

        const ReferenceElement& simplex () const
        {
          return reference_elements_[ Dune::Impl::SimplexTopology< dim >::type::id ];
        }

        const ReferenceElement& cube () const
        {
          return reference_elements_[ Dune::Impl::CubeTopology< dim >::type::id ];
        }

        const ReferenceElement& pyramid () const
        {
          return reference_elements_[ Dune::Impl::PyramidTopology< dim >::type::id ];
        }

        const ReferenceElement& prism () const
        {
          return reference_elements_[ Dune::Impl::PrismTopology< dim >::type::id ];
        }

        const_iterator begin () const
        {
          return reference_elements_.data();
        }

        const_iterator end () const
        {
          return reference_elements_.data() + numTopologies;
        }

      private:

        std::array<Implementation,numTopologies> implementations_;
        std::array<ReferenceElement,numTopologies> reference_elements_;

      };


    } // namespace Impl


#endif // DOXYGEN


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
    template< class ctype_, int dim >
    struct ReferenceElements
    {

      //! The coordinate field type of the contained reference elements.
      using ctype = ctype_;

      //! The coordinate field type of the contained reference elements.
      using CoordinateField = ctype;

      //! The dimension of the contained reference elements.
      static constexpr int dimension = dim;

    private:

      using Container = Impl::ReferenceElementContainer< ctype, dim >;

    public:

      //! The reference element type.
      using ReferenceElement = typename Container::ReferenceElement;

      //! Iterator over available reference elements.
      using Iterator         = typename Container::const_iterator;

      //! Iterator over available reference elements.
      using iterator         = Iterator;

      //! get general reference elements
      static const ReferenceElement&
      general ( const GeometryType& type )
      {
        return container() ( type );
      }

      //! get simplex reference elements
      static const ReferenceElement& simplex ()
      {
        return container().simplex();
      }

      //! get hypercube reference elements
      static const ReferenceElement& cube ()
      {
        return container().cube();
      }

      static Iterator begin ()
      {
        return container().begin();
      }

      static Iterator end ()
      {
        return container().end();
      }

    private:

      DUNE_EXPORT static const Container& container ()
      {
        static Container container;
        return container;
      }
    };

  } // namespace Geo

  using Geo::ReferenceElements;

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
