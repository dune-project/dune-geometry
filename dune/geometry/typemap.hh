// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_TYPEMAP_HH
#define DUNE_GEOMETRY_TYPEMAP_HH

#include <dune/common/array.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

namespace Dune
{

  // GeometryTypeMap
  // ---------------

  /** \brief associative container assigning values to each GeometryType
   *
   *  The GeometryTypeMap is an associative container similar to
   *  std::map< GeometryType, T >.
   *
   *  \tparam  mindim   minimum dimension to associate values to
   *  \tparam  maxdim   maximum dimension to associate values to
   *  \tparam  T        type of values
   *  \tparam  incnone  associate data to type none (defaults to true)
   *
   *  \note T must be default constructable and copyable
   */

  template< int mindim, int maxdim, class T, bool incnone = true >
  class GeometryTypeMap
  {
    typedef GeometryTypeMap< mindim, maxdim, T, incnone > This;

    typedef GeometryTypeIndex< mindim, maxdim, incnone > Index;
    typedef array< T, Index::size > Container;

  public:
    /** \brief type of values */
    typedef typename Container::value_type value_type;

    /** \brief type modeling a const reference (usually const values &) */
    typedef typename Container::const_reference const_reference;
    /** \brief type modeling a reference (usually values &) */
    typedef typename Container::reference reference;

    /** \brief type of iterator over constant list */
    typedef typename Container::const_iterator const_iterator;
    /** \brief type of iterator over mutable list */
    typedef typename Container::iterator iterator;
    /** \brief type of reverse iterator over constant list */
    typedef typename Container::const_reverse_iterator const_reverse_iterator;
    /** \brief type of reverse iterator over mutable list */
    typedef typename Container::reverse_iterator reverse_iterator;

    GeometryTypeMap () {}

    GeometryTypeMap ( const value_type &value ) { fill( value ); }

    /** \brief random access to data associated to a geometry type
     *
     *  \param[in]  gt  geometry type
     */
    const_reference operator[] ( const GeometryType &gt ) const { return container_[ Index::index( gt ) ]; }
    /** \brief random access to data associated to a geometry type
     *
     *  \param[in]  gt  geometry type
     */
    reference operator[] ( const GeometryType &gt ) { return container_[ Index::index( gt ) ]; }

    /** \brief obtain iterator pointing to the first element in the container */
    const_iterator begin () const { return container_.begin(); }
    /** \brief obtain iterator pointing to the first element in the container */
    iterator begin () { return container_.begin(); }
    /** \brief obtain iterator pointing behind the last element in the container */
    const_iterator end () const { return container_.end(); }
    /** \brief obtain iterator pointing behind the last element in the container */
    iterator end () { return container_.end(); }

    /** \brief obtain reverse iterator pointing behind the last element of the container */
    const_reverse_iterator rbegin () const { return container_.rbegin(); }
    /** \brief obtain reverse iterator pointing behind the last element of the container */
    reverse_iterator rbegin () { return container_.rbegin(); }
    /** \brief obtain reverse iterator pointing to the first element of the container */
    const_reverse_iterator rend () const { return container_.rend(); }
    /** \brief obtain reverse iterator pointing to the first element of the container */
    reverse_iterator rend () { return container_.rend(); }

    /** \brief assign the same value to all geometry types
     *
     *  \param[in]  value  value to fill the container with
     */
    void fill ( const value_type &value ) { container_.fill( value ); }
    /** \brief swap the content of two containers */
    void swap ( const This &other ) { container_.swap( other.container_ ); }

  private:
    Container container_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_TYPEMAP_HH
