#ifndef DUNE_GEOMETRY_REFERENCEELEMENT_HH
#define DUNE_GEOMETRY_REFERENCEELEMENT_HH

#include <dune/common/deprecated.hh>

#include <dune/geometry/type.hh>

namespace Dune {
  namespace Geo {

    namespace Impl {

      // forward declaration for friend declaration
      template<typename ctype, int dim>
      class ReferenceElementContainer;

    }

    // forward declaration for constructing default reference element type
    template<typename ctype, int dim>
    class ReferenceElementImplementation;

    // forward declaration for deprecation check reference element type
    template<typename ctype, int dim>
    class DeprecatedReferenceElement;

    // forward declaration for backwards compatibility conversion
    template<typename ctype, int dim>
    struct ReferenceElements;

    // ReferenceElement
    // ----------------

    /** \class ReferenceElement
     *  \ingroup GeometryReferenceElements
     *  \brief This class provides access to geometric and topological
     *  properties of a reference element.
     *
     *  This includes its type,
     *  the number of subentities, the volume, and a method for checking
     *  if a point is contained in the reference element.
     *  The embedding of each subentity into the reference element is also
     *  provided.
     *
     *  This class has value semantics, i.e. you can (and should) pass it
     *  around by value and not by reference and store a copy of it.
     *
     *  Instances of this object for a given geometry type can be retrieved
     *  from the ReferenceElements class.
     *
     */
    template<typename Implementation>
    class ReferenceElement
    {

    public:

      /** \brief Collection of types depending on the codimension */
      template<int codim>
      using Codim = typename Implementation::template Codim<codim>;

      //! The coordinate field type.
      using ctype = typename Implementation::ctype;

      //! The coordinate field type.
      using CoordinateField = ctype;

      //! The coordinate type.
      using Coordinate = typename Implementation::Coordinate;

      //! The dimension of the reference element.
      static constexpr int dimension = Implementation::dimension;


      /** \brief number of subentities of codimension c
       *
       *  \param[in]  c  codimension whose size is desired
       */
      int size(int c) const
      {
        return _impl->size(c);
      }


      /** \brief number of subentities of codimension cc of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the number of subentities
       *  of codimension cc of the current reference element, that are also
       *  a subentity of E.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  cc  codimension whose size is desired (c <= cc <= dim)
       */
      int size(int i, int c, int cc) const
      {
        return _impl->size(i,c,cc);
      }


      /** \brief obtain number of ii-th subentity with codim cc of (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. And denote by S the ii-th subentity of codimension
       *  (cc-c) of E. Then, S is a also a subentity of codimension c of the current
       *  reference element. This method returns the number of S with respect
       *  to the current reference element.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  ii  number of subentity S (with respect to E)
       *  \param[in]  cc  codimension of subentity S (c <= cc <= dim)
       */
      int subEntity(int i, int c, int ii, int cc) const
      {
        return _impl->subEntity(i,c,ii,cc);
      }


      /** \brief obtain the type of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the GeometryType of E.
       *
       *  \param[in]  i      number of subentity E (0 <= i < size( c ))
       *  \param[in]  c      codimension of subentity E
       */
      GeometryType type(int i, int c) const
      {
        return _impl->type(i,c);
      }


      /** \brief obtain the type of this reference element */
      GeometryType type() const
      {
        return _impl->type();
      }


      /** \brief position of the barycenter of entity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the coordinates of
       *  the center of gravity of E within the current reference element.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       */
      Coordinate position(int i, int c) const
      {
        return _impl->position(i,c);
      }


      /** \brief check if a coordinate is in the reference element
       *
       *  This method returns true if the given local coordinate is within this
       *  reference element.
       *
       *  \param[in]  local  coordinates of the point
       */
      bool checkInside(const Coordinate& local) const
      {
        return _impl->checkInside(local);
      }


      /** \brief obtain the embedding of subentity (i,codim) into the reference
       *         element
       *
       *  Denote by E the i-th subentity of codimension codim of the current
       *  reference element. This method returns a \ref Dune::AffineGeometry
       *  that maps the reference element of E into the current reference element.
       *
       *  \tparam     codim  codimension of subentity E
       *
       *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
       */
      template<int codim>
      typename Codim<codim>::Geometry geometry(int i) const
      {
        return _impl->template geometry<codim>(i);
      }


      /** \brief obtain the volume of the reference element */
      CoordinateField volume() const
      {
        return _impl->volume();
      }


      /** \brief obtain the integration outer normal of the reference element
       *
       *  The integration outer normal is the outer normal whose length coincides
       *  with the face's integration element.
       *
       *  \param[in]  face  index of the face, whose normal is desired
       */
      Coordinate integrationOuterNormal(int face) const
      {
        return _impl->integrationOuterNormal(face);
      }


      /** \brief Constructs an empty reference element.
       *
       * This constructor creates an empty (invalid) reference element. This element may not be
       * used in any way except for assigning other reference elements to it. After
       * assigning a valid reference element (obtained from ReferenceElements), it may
       * be used without restrictions.
       */
      ReferenceElement()
        : _impl(nullptr)
      {}

      /** \brief Returns a reference to the internal implementation object.
       *
       * \warning This method may only be called on valid reference elements.
       * \warning This method exposes undocumented internals that may change without notice!
       */
      const Implementation& impl() const
      {
        return *_impl;
      }

#ifndef DOXYGEN

      DUNE_DEPRECATED_MSG("Capturing reference elements by reference is deprecated in DUNE 2.6. Please store a copy instead.")
      operator const DeprecatedReferenceElement<ctype,dimension>&() const
      {
        return ReferenceElements<ctype,dimension>::deprecated(*this);
      }

#endif // DOXYGEN

    private:

      // The implementation must be a friend to construct a wrapper around itself.
      friend Implementation;

      // The reference container is a friend to be able to call setImplementation.
      friend class Impl::ReferenceElementContainer<ctype,dimension>;

      // Constructor for wrapping an implementation reference (required internally by the default implementation)
      ReferenceElement(const Implementation& impl)
        : _impl(&impl)
      {}

      void setImplementation(const Implementation& impl)
      {
        _impl = &impl;
      }

      const Implementation* _impl;

    };

  }

  // make type signature of wrapper compatible with old implementation
  template<typename ctype, int dim>
  using ReferenceElement = Geo::DeprecatedReferenceElement<ctype,dim>;


}


#endif // DUNE_GEOMETRY_REFERENCEELEMENT_HH
