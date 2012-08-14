// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_TYPEMAP_HH
#define DUNE_GEOMETRY_TYPEMAP_HH

#include <dune/common/container/array.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

namespace Dune
{

  // GeometryTypeMap
  // ---------------

  template< int mindim, int maxdim, class T, bool incnone = true >
  class GeometryTypeMap
  {
    typedef GeometryTypeMap< mindim, maxdim, T, incnone > This;

    typedef GeometryTypeIndex< mindim, maxdim, incnone > Index;
    typedef array< T, Index::size > Container;

  public:
    typedef typename Container::value_type value_type;

    typedef typename Container::const_reference const_reference;
    typedef typename Container::reference reference;

    typedef typename Container::const_iterator const_iterator;
    typedef typename Container::iterator iterator;
    typedef typename Container::const_reverse_iterator const_reverse_iterator;
    typedef typename Container::reverse_iterator reverse_iterator;

    const_reference operator[] ( const GeometryType &gt ) const { return container_[ Index::index( gt ) ]; }
    reference operator[] ( const GeometryType &gt ) { return container_[ Index::index( gt ) ]; }

    const_iterator begin () const { return container_.begin(); }
    iterator begin () { return container_.begin(); }
    const_iterator end () const { return container_.end(); }
    iterator end () { return container_.end(); }

    const_reverse_iterator rbegin () const { return container_.rbegin(); }
    reverse_iterator rbegin () { return container_.rbegin(); }
    const_reverse_iterator rend () const { return container_.rend(); }
    reverse_iterator rend () { return container_.rend(); }

    void fill ( const value_type &value ) { container_.fill( value ); }
    void swap ( const This &other ) { container_.swap( other.container_ ); }

  private:
    Container container_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_TYPEMAP_HH
