// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_PYTHON_GEOMETRY_MULTILINEARGEOMETRY_HH
#define DUNE_PYTHON_GEOMETRY_MULTILINEARGEOMETRY_HH

#include <type_traits>

#include <dune/common/visibility.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/python/common/typeregistry.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {
    // registerMemberFunctions
    // -----------------------

    template<typename MGeometry>
    struct RegisterMemberFunctions
    {
      RegisterMemberFunctions() {}
      void operator()(pybind11::class_<MGeometry>& cls)
      {
        cls.def("toLocal"                  , &MGeometry::local);
        cls.def("jacobianInverseTransposed", &MGeometry::jacobianInverseTransposed);
      }
    };



    // doNothing
    // ---------

    template<typename MGeometry>
    struct DoNothing
    {
      DoNothing() {}
      void operator()(pybind11::class_<MGeometry>&) {}
    };



    // registerMultiLinearGeometryType
    // -------------------------------

    template<class ctype, int dim, int coorddim>
    inline auto registerMultiLinearGeometryType(pybind11::module scope)
    {
      using pybind11::operator""_a;

      typedef MultiLinearGeometry<ctype, dim, coorddim> MGeometry;

      auto entry = insertClass< MGeometry >( scope, "MultiLinearGeometry_"+std::to_string(dim)+"_"+std::to_string(coorddim),
          GenerateTypeName("MultiLinearGeometry",MetaType<ctype>(),dim,coorddim),
          IncludeFiles{"dune/geometry/multilineargeometry.hh"}
          );
      auto cls = entry.first;

      if (!entry.second)
      {
        cls.def( pybind11::init( [] ( Dune::GeometryType type, pybind11::list corners ) {
            const std::size_t numCorners = corners.size();
            std::vector< FieldVector< double, coorddim > > copyCorners( numCorners );
            for( std::size_t i = 0; i < numCorners; ++i )
              copyCorners[ i ] = corners[ i ].template cast< FieldVector< double, coorddim > >();
            return new MGeometry( type, copyCorners );
          } ), "type"_a, "corners"_a );

        // toLocal and jacobianInverseTransposed call
        // MatrixHelper::template (xT)rightInvA<m, n> where n has to be >= m (static assert)
        std::conditional_t<(coorddim >= dim), RegisterMemberFunctions<MGeometry>, DoNothing<MGeometry> >{}(cls);

        cls.def_property_readonly("affine" , [](const MGeometry& geom) { return geom.affine(); });
        cls.def_property_readonly("type"   , &MGeometry::type);
        cls.def_property_readonly("corners", &MGeometry::corners);
        cls.def_property_readonly("center" , &MGeometry::center);
        cls.def_property_readonly("volume" , &MGeometry::volume);
        cls.def("corner"                   , &MGeometry::corner);
        cls.def("integrationElement"       , &MGeometry::integrationElement);

        cls.def("jacobianTransposed",
            [](const MGeometry& geom, const typename MGeometry::LocalCoordinate& local) {
              return geom.jacobianTransposed(local);
            });

        cls.def("toGlobal",
            [](const MGeometry& geom, const typename MGeometry::LocalCoordinate& local) {
                return geom.global(local);
              });
      }
#if 0
      else
      {
        scope.def( detail::nameMultiLinearGeometry< dim, coorddim >, [] ( Dune::GeometryType gt, pybind11::list corners ) {
            std::vector<FieldVector<double, 1>> cornerValues(corners.size());
            for (unsigned i = 0; i < corners.size(); ++i)
              cornerValues[i] = corners[i].cast<double>();
            return MGeometry(gt, cornerValues);
          } );
      }
#endif
      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // ifndef DUNE_PYTHON_GEOMETRY_MULTILINEARGEOMETRY_HH
