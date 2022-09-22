// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTS_HH

#include <cassert>

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include <array>

#include <dune/common/typetraits.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/referenceelementimplementation.hh>

namespace Dune
{

  namespace Geo
  {

#ifndef DOXYGEN


    template<typename ctype, int dim>
    class DeprecatedReferenceElement
      : public ReferenceElement<ReferenceElementImplementation<ctype,dim>>
    {

    protected:

      DeprecatedReferenceElement() = default;

    public:

      DeprecatedReferenceElement(const DeprecatedReferenceElement&) = delete;
      DeprecatedReferenceElement& operator=(const DeprecatedReferenceElement&) = delete;

      DeprecatedReferenceElement(const ReferenceElementImplementation<ctype,dim>& impl)
        : ReferenceElement<ReferenceElementImplementation<ctype,dim>>(impl)
      {}

    };


    template<typename ctype, int dim>
    class ConstructibleDeprecatedReferenceElement
      : public DeprecatedReferenceElement<ctype,dim>
    {
    public:
      ConstructibleDeprecatedReferenceElement() = default;
    };


    namespace Impl
    {

      // ReferenceElementContainer
      // -------------------------

      template< class ctype, int dim >
      class ReferenceElementContainer
      {
        static const unsigned int numTopologies = dim >= 0 ? (1u << dim) : 0;

        using Implementation   = ReferenceElementImplementation< ctype, dim >;
        using ConstructibleDeprecatedReferenceElement = Dune::Geo::ConstructibleDeprecatedReferenceElement<ctype,dim>;

      public:

        using DeprecatedReferenceElement = Dune::Geo::DeprecatedReferenceElement<ctype,dim>;

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
          return reference_elements_[ Dune::GeometryTypes::simplex(dim).id() ];
        }

        const ReferenceElement& cube () const
        {
          return reference_elements_[ Dune::GeometryTypes::cube(dim).id() ];
        }

        const ReferenceElement& pyramid () const
        {
          return reference_elements_[ Dune::GeometryTypes::pyramid.id() ];
        }

        const ReferenceElement& prism () const
        {
          return reference_elements_[ Dune::GeometryTypes::prism.id() ];
        }

        const_iterator begin () const
        {
          return reference_elements_.data();
        }

        const_iterator end () const
        {
          return reference_elements_.data() + numTopologies;
        }

        // here, we make sure to actually return a const reference to something
        // that is guaranteed not to become invalid, as otherwise, we might run
        // straight into debugging hell when a user binds the return value to a
        // const ref and the temporary goes out of scope.
        const DeprecatedReferenceElement& deprecated(const ReferenceElement& v) const
        {
          return reference_elements_[v.impl().type(0,0).id()];
        }

      private:

        std::array<Implementation,numTopologies> implementations_;
        std::array<ConstructibleDeprecatedReferenceElement,numTopologies> reference_elements_;

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

#ifndef DOXYGEN
      static const typename Container::DeprecatedReferenceElement&
      deprecated(const ReferenceElement& v)
      {
        return container().deprecated(v);
      }
#endif // DOXYGEN

    private:

      DUNE_EXPORT static const Container& container ()
      {
        static Container container;
        return container;
      }
    };

  } // namespace Geo

  //! Make the ReferenceElements container available in the old location.
  using Geo::ReferenceElements;


#ifdef DOXYGEN

  //! Returns a reference element for the objects t....
  /**
   * \ingroup GeometryReferenceElements
   * The freestanding function referenceElement is a generic entry point
   * for getting reference elements for arbitrary objects that support the
   * operation. As it relies on argument-dependent lookup, the function should
   * be called without any qualifying namespace. Note, however, that the versions
   * of referenceElement() for a dimension and GeometryType with explicit template
   * arguments cannot be found by ADL, so you have to explicitly pull them in or
   * qualify the call using `Dune::`:
   *
     \code
     {
       // option 1: using
       using Dune::referenceElement;
       auto ref_el = referenceElement<double,dim>(geometry_type);
     }
     {
       // option 2: explicitly put in Dune::
       auto ref_el = Dune::referenceElement<double,dim>(geometry_type);
     }
     {
       // option 3: use version without explicit template arguments
       auto ref_el = referenceElement(double(),geometry_type,Dune::Dim<dim>());
     }
     \endcode
   *
   * The returned object is guaranteed to have value semantics, so you can copy it
   * around and store it by value.
   *
   * The grid geometries in dune-grid support this function, and thus most people
   * will use this function as
   *
     \code
     for (const auto& cell : elements(grid_view))
       {
         auto geo = cell.geometry();
         auto ref_el = referenceElement(geo);
         // do some work...
       }
     \endcode
   *
   * This does of course also work for entities of other codimensions.
   */
  template<typename... T>
  unspecified-value-type referenceElement(T&&... t);

#endif


  //! Returns a reference element of dimension dim for the given geometry type and coordinate field type.
  /**
   * This function allows you to obtain a reference element for a given coordinate type, dimension and
   * GeometryType:
   *
     \code
     auto gt = ...;
     using Dune::referenceElement;
     auto ref_el = referenceElement<ctype,dim>(gt);
     auto ref_el = referenceElement<ctype>(gt,Dune::Dim<dim>());
     \endcode
   *
   * \ingroup GeometryReferenceElements
   */
  template<typename T, int dim>
  auto referenceElement(const Dune::GeometryType& gt, Dune::Dim<dim> = {})
  {
    return ReferenceElements<T,dim>::general(gt);
  }


  //! Returns a reference element of dimension dim for the given geometry type and coordinate field type.
  /**
   * This function allows you to obtain a reference element for a given coordinate type, dimension and
   * GeometryType:
   *
     \code
     auto gt = ...;
     using Dune::referenceElement;
     auto ref_el = referenceElement(ctype(),gt,Dune::Dim<dim>());
     \endcode
   *
   * \ingroup GeometryReferenceElements
   */
  template<typename T, int dim, std::enable_if_t<IsNumber<std::decay_t<T>>::value, int> = 0>
  auto referenceElement(const T&, const Dune::GeometryType& gt, Dune::Dim<dim>)
  {
    return ReferenceElements<T,dim>::general(gt);
  }


#ifndef DOXYGEN

  // helpers for the ReferenceElement<> meta function

  namespace Impl {

    // Evaluates to the correct reference element iff <T...> matches the pattern <number_type,Dim<int>>
    // otherwise, it's ill-formed. Should be used with detected_or and friends.

    template<typename... T>
    struct DefaultReferenceElementExtractor;

    template<typename T, typename std::enable_if<IsNumber<T>::value,int>::type dim>
    struct DefaultReferenceElementExtractor<T,Dim<dim>>
    {
      using type = typename Dune::Geo::ReferenceElements<T,dim>::ReferenceElement;
    };

    template<typename... T>
    using DefaultReferenceElement = typename DefaultReferenceElementExtractor<T...>::type;

  }

  // looks up the type of a reference element by trying to instantiate the correct overload
  // of referenceElement() for the given arguments. This will fail if there is no valid
  // overload and should be used with detected_or or some other utility that places the
  // instantiation in SFINAE context.
  //
  // this is placed directly in namespace Dune to avoid any weird surprises

  template<typename... T>
  using LookupReferenceElement = decltype(referenceElement(std::declval<T>()...));

#endif // DOXYGEN

  namespace Transitional {

#ifndef DOXYGEN

    // this abomination checks whether the template signature matches the special case
    // ReferenceElement<number_type,Dune::Dim<int>> and otherwise defers the type lookup
    // to a decltype on a call to referenceElement(std::declval<T>())

    template<typename... T>
    using ReferenceElement = Std::detected_or_t<
      Std::detected_t<LookupReferenceElement,T...>,
      Impl::DefaultReferenceElement,
      T...
      >;

#else // DOXYGEN

    //! Returns the type of reference element for the argument types T...
    /**
     * This type alias can be used to get the type of a reference element if you
     * want to store the element as a class member or need to access nested type
     * information. Normally, it will return the type of reference element that a
     * call to `referenceElement(T...)` would return, for example for a geometry:
     *
       \code
       Cell::Geometry geo = cell.geometry();
       Dune::ReferenceElement<Cell::Geometry> ref_el = referenceElement(geo);
       \endcode
     *
     * There is also a special shorthand signature for the default reference elements:
     *
       \code
       constexpr int dim = ...;
       auto geometry_type = ...;
       Dune::ReferenceElement<double,Dune::Dim<dim>> = referenceElement<double,dim>(geometry_type);
       \endcode
     *
     * \ingroup GeometryReferenceElements
     * \sa Dune::ReferenceElement
     */
    template<typename... T>
    using ReferenceElement = unspecified-type;

#endif //DOXYGEN

  }

#ifndef DOXYGEN

  namespace Impl {

    // helpers for triggering a deprecation warning for occurrences of the old
    // ReferenceElement syntax (Dune::ReferenceElement<T,int>)

    // this looks a little weird, encoding the type in the return type of a nested function,
    // but it was the only reliable way to trigger the warning iff the compiler encounters
    // an invalid use. Other solutions either miss some (or all) uses or trigger false alarms.

    template<typename T>
    struct ValidReferenceElementTypeSignature
    {
      Transitional::ReferenceElement<T> check();
    };

    template<typename T>
    struct DeprecatedReferenceElementTypeSignature
    {
      [[deprecated("Dune::ReferenceElement<T,int> is deprecated, please use Dune::ReferenceElement<Geometry> (if you have a geometry), Dune::ReferenceElements<T,int>::ReferenceElement or Dune::Transitional::ReferenceElement<T,Dune::Dim<int>> instead. After Dune 2.6, you will be able to use Dune::ReferenceElement<T,Dune::Dim<int>>.")]] T check();
    };

  } // namespace Impl

  // This just makes sure the user doesn't use invalid syntax (by checking against the -1 default for
  // the dimension, and then either hands back the old-style type along with a deprecation warning or
  // forwards to the transitional version
  template<typename T, int dim = -1>
  using ReferenceElement = decltype(
    std::declval<
      typename std::conditional<
        dim == -1,
        Impl::ValidReferenceElementTypeSignature<T>,
        Impl::DeprecatedReferenceElementTypeSignature<Geo::DeprecatedReferenceElement<T,dim>>
        >::type
      >().check());

#else // DOXYGEN

  //! Returns the type of reference element for the argument type T.
  /**
   * This type alias can be used to get the type of a reference element if you
   * want to store the element as a class member or need to access nested type
   * information. Normally, it will return the type of reference element that a
   * call to `referenceElement(T)` would return, for example for a geometry:
   *
     \code
     Cell::Geometry geo = cell.geometry();
     Dune::ReferenceElement<Cell::Geometry> ref_el = referenceElement(geo);
     \endcode
   *
   * In the long run, we want to support multiple types T here, but for now that
   * does not work due to backwards compatibility reasons. You can also still obtain
   * the type of a standard reference element using Dune::ReferenceElement<ctype,dim>,
   * but this is deprecated in DUNE 2.6.
   * If you need support for ReferenceElement with multiple type arguments, you can use
   * Dune::Transitional::ReferenceElement for now, which supports multiple types, but not
   * the backwards compatibility mode and will become the default after the release of DUNE 2.6.
   * There is also a special shorthand signature for the default reference elements:
   *
   * \deprecated Using the syntax Dune::ReferenceElement<ctype,dim> is deprecated in DUNE 2.6.
   *             You have the following alternatives:
   *             - Most of the time, you will want to use Dune::ReferenceElement<Geometry> because
   *               you already have a geometry type available.
   *             - Dune::ReferenceElements<ctype,dim>::ReferenceElement.
   *             - Dune::Transitional::ReferenceElement<ctype,Dune::Dim<dim>>. This will become
   *               available as Dune::ReferenceElement<ctype,Dune::Dim<dim>> after the release of
   *               DUNE 2.6.
   *
   * \ingroup GeometryReferenceElements
   * \sa Dune::Transitional::ReferenceElement
   */
  template<typename T, int dim>
  using ReferenceElement = unspecified-type;

#endif // DOXYGEN



} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTS_HH
