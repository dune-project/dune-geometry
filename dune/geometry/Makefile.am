geometryincludedir = $(includedir)/dune/geometry
geometryinclude_HEADERS =			\
	affinegeometry.hh			\
	axisalignedcubegeometry.hh		\
	generalvertexorder.hh			\
	genericreferenceelements.hh		\
	mockgeometry.hh				\
	multilineargeometry.hh			\
	quadraturerules.hh			\
	referenceelements.hh			\
	refinement.hh				\
	topologyfactory.hh			\
	type.hh					\
	typeindex.hh				\
	virtualrefinement.hh			\
	virtualrefinement.cc

SUBDIRS = quadraturerules genericgeometry refinement test

noinst_LTLIBRARIES = libgeometry.la
libgeometry_la_SOURCES =
libgeometry_la_LIBADD = \
	quadraturerules/libquadraturerules.la \
	genericgeometry/libgenericgeometry.la \
	$(DUNE_LIBS)

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
