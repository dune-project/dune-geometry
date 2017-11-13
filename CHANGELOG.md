# Master (will become release 2.7)

- The class `GeometryType` has been cleaned up in major way:

  - The class and most of its methods are now `constexpr`.

  - There are new singletons and factory functions in the namespace `Dune::GeometryTypes`. These
    are now the official way to obtain a `GeometryType`.

  - `GeometryType::BasicType` and the assorted constructor have been deprecated and will be removed
    after the release of DUNE 2.6.

  - The assorted member functions `GeometryType::make...()` have been deprecated and will be removed
    after the release of DUNE 2.6.
