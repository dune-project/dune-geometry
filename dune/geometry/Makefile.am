noinst_LTLIBRARIES = libgridcommon.la
# construct only from sub-libraries
libgridcommon_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libgridcommon_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libgridcommon_la_LIBADD = quadraturerules/libquadraturerules.la \
  $(DUNE_LIBS)

geometryincludedir = $(includedir)/dune/geometry
geometryinclude_HEADERS = genericreferenceelements.hh geometry.hh \
	quadraturerules.hh topologyfactory.hh type.hh

SUBDIRS = quadraturerules genericgeometry

include $(top_srcdir)/am/global-rules
