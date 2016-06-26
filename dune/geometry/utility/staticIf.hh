// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_STATICIF_HH
#define DUNE_LOCALFUNCTIONS_COMMON_STATICIF_HH

#include <type_traits>

namespace Dune
{
  // implementation of staticIf
  namespace Impl
  {
    template <bool condition>
    struct StaticIf;

    template <>
    struct StaticIf<true>
    {
      // applies \param args... to the \param then_ branch only
      template <class Then, class Else, class... Args>
      static auto eval(Then&& then_, Else&&, Args&&... args)
      {
        return then_(std::forward<Args>(args)...);
      }
    };

    template <>
    struct StaticIf<false>
    {
      // applies \param args... to the \param else_ branch only
      template <class Then, class Else, class... Args>
      static auto eval(Then&&, Else&& else_, Args&&... args)
      {
        return else_(std::forward<Args>(args)...);
      }
    };

  } // end namespace Impl



  /**
   * \brief Static if
   *
   * \ingroup Utility
   *
   * Calls the branch \param then_ or \param else_ depending on a compile-time
   * condition. The branches can be callable with a List of arguments of type
   * \param args... An example of usage is:
   * ```
   * staticIf< Dune::Concept::localBasisHasPartial<Imp>() >
   * (
   *    [&](auto&& impl) { impl.partial(order, in, out); },    // then_
   *    [] (auto&& impl) { DUNE_THROW(NotImplemented, "!"); }, // else_
   *    impl_                                                  // args...
   * );
   * ```
   * This calls the partial() method of the argument `impl` only, if the
   * corresponding class has such a method and throws and exception otherwise.
   * Thus, a compile-time error can be transfered into a runtime error.
   */
  template <bool condition, class Then, class Else, class... Args>
  auto staticIf(Then&& then_, Else&& else_, Args&&... args)
  {
    return Impl::StaticIf<condition>::eval(std::forward<Then>(then_),
                                           std::forward<Else>(else_),
                                           std::forward<Args>(args)...);
  }

  /// Same as \ref staticIf, but ignores the else_ branch.
  template <bool condition, class Then>
  void staticIf(Then&& then_)
  {
    Impl::StaticIf<condition>::eval(std::forward<Then>(then_), [](){ /* do nothing */ } );
  }

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_COMMON_STATICIF_HH
