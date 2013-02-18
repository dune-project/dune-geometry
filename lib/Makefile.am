# $Id: $
dunemodulelibdir=$(libdir)

# the dune-cgeometry library
lib_LTLIBRARIES = libdunegeometry.la

# this is just a renaming of libgeometry
libdunegeometry_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libdunegeometry_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libdunegeometry_la_LIBADD = ../dune/geometry/libgeometry.la

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
