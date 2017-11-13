# Release 2.6

- The enum `GeometryType::BasicType` is deprecated, and will be removed after Dune 2.6.

- `VirtualRefinement` and `Refinement` now support arbitrary refinements, not
  just powers of two.  Whereever you where passing a parameter `int levels`
  (now deprecated), you should now pass a parameter `RefinementIntervals
  intervals`.  There are convenience functions `refinementIntervals(int
  intervals)` and `refinementLevels(int levels)` to construct parameters of
  type `RefinementIntervals`.

- The class `GeometryType` has been cleaned up in major way:

  - The class and most of its methods are now `constexpr`.

  - There are new singletons and factory functions in the namespace `Dune::GeometryTypes`. These
    are now the official way to obtain a `GeometryType`.

  - `GeometryType::BasicType` and the assorted constructor have been deprecated and will be removed
    after the release of DUNE 2.6.

  - The assorted member functions `GeometryType::make...()` have been deprecated and will be removed
    after the release of DUNE 2.6.

- The reference element interface has had a substantial overhaul that can break backwards
  compatibility in some corner cases.

  - `ReferenceElement` has value semantics now: You should store instances by value and can freely
    copy them around. Doing so is not more expensive than storing a const reference.

  - As a consequence of value semantics, `ReferenceElement` is default constructible now. A default
    constructed `ReferenceElement` may only be assigned another `ReferenceElement`; all other
    oeprations cause undefined behavior. Moreover, instances are now comparable and hashable to
    allow storing them in maps.

  - We have added code that tries to warn you if you are still storing a `ReferenceElement` by const
    reference; please update all those occurences.

  - The meaning of `Dune::ReferenceElement` has changed. It is not a type anymore, but an alias
    template that looks up the correct implementation for the given template arguments. For now,
    there is only a single implementation, but we expect people to come up with additional
    implementations in the future. For this reason, the syntax `Dune::ReferenceElement<ctype,dim>`
    is deprecated and will cause compilation failures in the future. If you still need access to
    that type, use `typename Dune::ReferenceElements<ctype,dim>::ReferenceElement` instead.

  - You can now directly obtain a reference element for a given geometry using the free function
    `referenceElement(geometry)`. This function should be called without any namespace qualifiers to
    enable ADL and you should normally capture the return value of the function using `auto`, but if
    you need to explicitely access the type, this is also available as
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
