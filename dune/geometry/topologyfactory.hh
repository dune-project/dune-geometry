// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TOPOLOGYFACTORY_HH
#define DUNE_GEOMETRY_TOPOLOGYFACTORY_HH

#include <cassert>

#include <array>
#include <map>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

namespace Dune
{

  /**
   * @brief Provide a factory over the generic topologies
   *
   * This class can be used to dynamically create objects
   * statically bound by their generic topologies.
   * The method create() returns a pointer to an object depending
   * on the topology id and a key; the dimension corresponding
   * to the topology id is static and is provided by the
   * Traits class. A static method (taking the Topology as template
   * argument) is also provided.
   * The Traits class must provide the space dimension,
   * the types for the key (Key),
   * the objects returned (Object), and the underlying factory
   * (Factory). This class must have a template method
   * createObject taking a key and returning a pointer to
   * the newly create Object - for destruction call the release
   * method.
   **/
  template <class Traits>
  struct TopologyFactory
  {
    // extract types from Traits class
    static const unsigned int dimension = Traits::dimension;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;
    typedef typename Traits::Factory Factory;

    //! dynamically create objects
    static Object *create ( const Dune::GeometryType &gt, const Key &key )
    {
      return Impl::toGeometryTypeIdConstant<dimension>(gt, [&](auto id) {
        return create<decltype(id)::value>(key);
      });
    }
    //! statically create objects
    template< GeometryType::Id geometryId >
    static Object *create ( const Key &key )
    {
      return Factory::template createObject< geometryId >( key );
    }

    //! statically create objects
    template< class Topology >
    static Object *create ( const Key &key )
    {
      return Factory::template createObject< Topology >( key );
    }

    //! release the object returned by the create methods
    static void release( Object *object ) { delete object; }
  };



  /** @brief A wrapper for a TopologyFactory providing
   *         singleton storage. Same usage as TopologyFactory
   *         but with empty release method an internal storage.
   **/
  template <class Factory>
  struct TopologySingletonFactory
  {
    static const unsigned int dimension = Factory::dimension;
    typedef typename Factory::Key Key;
    typedef const typename Factory::Object Object;

    //! @copydoc TopologyFactory::create(const Dune::GeometryType &gt,const Key &key)
    static Object *create ( const Dune::GeometryType &gt, const Key &key )
    {
      assert( gt.id() < numTopologies );
      return instance().getObject( gt, key );
    }
    //! @copydoc TopologyFactory::create(const Key &key)
    template< GeometryType::Id geometryId >
    static auto create ( const Key &key )
      -> std::enable_if_t< static_cast<GeometryType>(geometryId).dim() == dimension, Object * >
    {
      return instance().template getObject< geometryId >( key );
    }

    //! @copydoc TopologyFactory::create(const Key &key)
    template< class Topology >
    static auto create ( const Key &key )
      -> std::enable_if_t< Topology::dimension == dimension, Object * >
    {
      return instance().template getObject< Topology >( key );
    }

    //! @copydoc TopologyFactory::release
    static void release ( Object *object )
    {}

  private:
    struct ObjectDeleter
    {
      void operator() ( Object *ptr ) const { Factory::release( ptr ); }
    };

    static TopologySingletonFactory &instance ()
    {
      static TopologySingletonFactory instance;
      return instance;
    }

    static const unsigned int numTopologies = (1 << dimension);
    typedef std::array< std::unique_ptr< Object, ObjectDeleter >, numTopologies > Array;
    typedef std::map< Key, Array > Storage;

    TopologySingletonFactory () = default;

    std::unique_ptr< Object, ObjectDeleter > &find ( const unsigned int topologyId, const Key &key )
    {
      return storage_[ key ][ topologyId ];
    }

    Object *getObject ( const Dune::GeometryType &gt, const Key &key )
    {
      auto &object = find( gt.id(), key );
      if( !object )
        object.reset( Factory::create( gt, key ) );
      return object.get();
    }

    template< GeometryType::Id geometryId >
    Object *getObject ( const Key &key )
    {
      static constexpr GeometryType geometry = geometryId;
      auto &object = find( geometry.id(), key );
      if( !object )
        object.reset( Factory::template create< geometry >( key ) );
      return object.get();
    }

    template< class Topology >
    Object *getObject ( const Key &key )
    {
      auto &object = find( Topology::id, key );
      if( !object )
        object.reset( Factory::template create< Topology >( key ) );
      return object.get();
    }

    Storage storage_;
  };

}

#endif // #ifndef DUNE_GEOMETRY_TOPOLOGYFACTORY_HH
