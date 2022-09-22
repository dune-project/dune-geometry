<!--
SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
-->

# Master (will become release 2.9)

- The `Geometry` interface was extended by methods `jacobian(local)` and `jacobianInverse(local)`
  and corresponding typedefs `Jacobian` and `JacobianInverse`. This is implemented by all geometry
  implementations provided by dune-geometry. But external implementations need to be adjusted
  to pass the interface check provided by `checkgeometry.hh`.

- The `Geometry::integrationElement` now needs to return the type `Volume`
  instead of `ctype`. Note that this may be different from `ctype` if the
  geometry supports typed dimensions. In such case, `ctype` is a length, and not
  appropriate for a volume quantity.

# Release 2.8

- Python bindings have been moved from the `dune-python` module which is now
  obsolete. To activate Python bindings the CMake flag
  `DUNE_ENABLE_PYTHONBINDINGS` needs to be turned on (default is off).
  Furthermore, flags for either shared library or position independent code
  needs to be used.

- The class `AxisAlignedCubeGeometry` has always had a constructor taking
  two arguments `FieldVector<ctype,coorddim> lower` and `FieldVector<ctype,coorddim> upper`.
  This constructor was always to be used in the case `dim==coorddim` only,
  but this was never enforced.  Starting with version 2.8, compilation
  fails with an error message if this constructor is used with `dim!=coorddim`.

- Two new sets of quadrature rules are provided: the left and right Gauss-Radau quadrature rules.
  These are optimal rules that include only one endpoint of the integration interval
  (either left or right) and integrate polynomials of order 2n - 2 exactly.

- GeometryType has four new methods: `isPrismatic()`, `isPrismatic(int step)` and `isConical()`,`isConical(int step)`.
  The versions with an argument return true if the corresponding construction was used in step 0 <= `step` <=dim-1.
  The other two assume a default argument of `dim-1` (the latest construction step).

- GeometryTypes has two new methods: `prismaticExtension(GeometryType gt)` and `conicalExtension(GeometryType gt)`.
  They return an extended GeometryType based on `gt` via the corresponding construction. For example:
  ```c++
  GeometryType gt = GeometryTypes::line;
  auto square = GeometryTypes::prismaticExtension(gt);
  auto triangle = GeometryTypes::conicalExtension(gt);
  ```

## Deprecations and removals

- Remove code needed to use reference elements by reference.

- Remove `GeometryType`'s deprecated member functions
  `GeometryType::make...()`.

- Remove deprecated constructor `GeometryType(unsigned dim)`.

- Remove deprecated `CompositeQuadratureRule(QuadratureRule, int)`. Use
  `CompositeQuadratureRule(QuadratureRule, Dune::refinement{Intervals|Levels}(int))`
  instead.

- Removed all structs from `Impl` dealing with the recursive topology construction: `TopologyFactory`, `TopologySingletonFactory`,
  `Point`, `Prism`, `Pyramid`, `IsSimplex`, `IsCube`, `SimplexTopology`, `CubeTopology`, `PyramidTopology`, `PrismTopology`, `IfTopology`.
  Deprecated the free function `Impl::isTopology`.
  Use the geometries provided by `GeometryType` and `GeometryTypes` instead.
  To simplify the transition you can include the header "dune/geometry/deprecated_topology.hh".

# Release 2.7

- The reference elements have a new method `subEntities`. The result of
  `referenceELement.subEntities(i,codim, c)` is an iterable range
  containing the indices of all codim-`c` subentities of the subentity
  `(i,codim)`, e.g., the vertices of an edge. The range also provides
  the methods `size()` and `contains()`.
- The methods `GeometryType(int)` and `GeometryType(unsigned int)` have been deprecated
  and will be removed after the release of dune-geometry 2.7.  Instead, please now use
  `GeometryTypes::cube(dim)` to construct one- or two-dimensional `GeometryType` objects.
- Geometry implementations now export a type `Volume` that is used for the return
  value of the `volume` methods.  So does the generic `ReferenceElement` implementation.
-   More efficient quadrature rules for simplices are available that
    need less quadrature points to achieve the same order.  For now these
    have to be explicitly requested:
    ```c++
    auto&& rule = Dune::QuadratureRules<...>::rule(..., Dune::QuadratureType::GaussJacobi_n_0);
    ```
    See [!127].

    [!127]: https://gitlab.dune-project.org/core/dune-geometry/merge_requests/127

# Release 2.6

- The enum `GeometryType::BasicType` is deprecated, and will be removed after Dune 2.6.

- `VirtualRefinement` and `Refinement` now support arbitrary refinements, not
  just powers of two.  Wherever you where passing a parameter `int levels`
  (now deprecated), you should now pass a parameter `RefinementIntervals
  intervals`.  There are convenience functions `refinementIntervals(int
  intervals)` and `refinementLevels(int levels)` to construct parameters of
  type `RefinementIntervals`.

    See core/dune-geometry!51

- The class `GeometryType` has been cleaned up in major way:

    See core/dune-geometry!64 and core/dune-geometry!55

  - The class and most of its methods are now `constexpr`.

  - There are new singletons and factory functions in the namespace `Dune::GeometryTypes`. These
    are now the official way to obtain a `GeometryType`.

  - The constructor taking a `GeometryType::BasicType` and a dimension has been deprecated and will be
    removed after the release of DUNE 2.6.

  - The assorted member functions `GeometryType::make...()` have been deprecated and will be removed
    after the release of DUNE 2.6.

- The reference element interface has had a substantial overhaul that can break backwards
  compatibility in some corner cases.

    See core/dune-geometry!52

  - `ReferenceElement` has value semantics now: You should store instances by value and can freely
    copy them around. Doing so is not more expensive than storing a const reference.

  - As a consequence of value semantics, `ReferenceElement` is default constructible now. A default
    constructed `ReferenceElement` may only be assigned another `ReferenceElement`; all other
    operations cause undefined behavior. Moreover, instances are now comparable and hashable to
    allow storing them in maps.

  - We have added code that tries to warn you if you are still storing a `ReferenceElement` by const
    reference; please update all those occurrences.

  - The meaning of `Dune::ReferenceElement` has changed. It is not a type anymore, but an alias
    template that looks up the correct implementation for the given template arguments. For now,
    there is only a single implementation, but we expect people to come up with additional
    implementations in the future. For this reason, the syntax `Dune::ReferenceElement<ctype,dim>`
    is deprecated and will cause compilation failures in the future. If you still need access to
    that type, use `typename Dune::ReferenceElements<ctype,dim>::ReferenceElement` instead.

  - You can now directly obtain a reference element for a given geometry using the free function
    `referenceElement(geometry)`. This function should be called without any namespace qualifiers to
    enable ADL and you should normally capture the return value of the function using `auto`, but if
    you need to explicitly access the type, this is also available as
    `Dune::ReferenceElement<Geometry>`.

      In short: If you can, use the following idiom to obtain a reference element for a geometry:
      ```c++
      auto ref_el = referenceElement(geometry);
      ```

      The change to the meaning of `Dune::ReferenceElement` can break compilation if you have function
      overloads that partially specialize on it, e.g.
      ```c++
      template<typename ctype, int dim>
      void f(const Dune::ReferenceElement<ctype,dim> ref_el)
      {}
      ```
      Normally, you can just simplify this to the following code that also shows how to extract the
      missing template parameters:
     ```c++
      template<typename RefEl>
      void f(const RefEl ref_el)
      {
        using ctype = typename RefEl::CoordinateField;
        constexpr auto dim = RefEl::dimension;
      }
      ```
