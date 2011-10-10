noinst_LTLIBRARIES = libgeometry.la
# construct only from sub-libraries
libgeometry_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libgeometry_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libgeometry_la_LIBADD = quadraturerules/libquadraturerules.la \
  $(DUNE_LIBS)

geometryincludedir = $(includedir)/dune/geometry
geometryinclude_HEADERS = genericreferenceelements.hh  mockgeometry.hh \
	quadraturerules.hh  topologyfactory.hh  type.hh  vertexorder.hh

SUBDIRS = quadraturerules genericgeometry

include $(top_srcdir)/am/global-rules
