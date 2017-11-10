# Release 2.6

- The enum `GeometryType::BasicType` is deprecated, and will be removed after Dune 2.6.

- `VirtualRefinement` and `Refinement` now support arbitrary refinements, not
  just powers of two.  Whereever you where passing a parameter `int levels`
  (now deprecated), you should now pass a parameter `RefinementIntervals
  intervals`.  There are convenience functions `refinementIntervals(int
  intervals)` and `refinementLevels(int levels)` to construct parameters of
  type `RefinementIntervals`.
