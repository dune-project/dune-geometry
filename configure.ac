# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common
# this implies checking for [dune-common]
DUNE_CHECK_ALL

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS)')
AC_SUBST([LDADD], '$(top_builddir)/lib/libdunegeometry.la $(DUNE_LIBS)')

AC_CONFIG_FILES([
  Makefile
  lib/Makefile
  cmake/Makefile
  cmake/pkg/Makefile
  doc/Makefile
  doc/appl/Makefile
  doc/appl/refelements/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  doc/refinement/Makefile
  dune/Makefile
  dune/geometry/Makefile
  dune/geometry/test/Makefile
  dune/geometry/genericgeometry/Makefile
  dune/geometry/genericgeometry/test/Makefile
  dune/geometry/quadraturerules/Makefile
  dune/geometry/refinement/Makefile
  m4/Makefile
  dune-geometry.pc
])
AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
