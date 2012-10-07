geometryincludedir = $(includedir)/dune/geometry
geometryinclude_HEADERS =			\
	affinegeometry.hh			\
        axisalignedcubegeometry.hh              \
	bilineargeometry.hh			\
	generalvertexorder.hh			\
	multilineargeometry.hh			\
	mockgeometry.hh				\
	quadraturerules.hh			\
	referenceelements.hh			\
	topologyfactory.hh			\
	trilineargeometry.hh			\
	type.hh					\
	typeindex.hh				\
	typemap.hh

SUBDIRS = quadraturerules genericgeometry test

noinst_LTLIBRARIES = libgeometry.la
libgeometry_la_SOURCES = type.cc
libgeometry_la_LIBADD = \
	quadraturerules/libquadraturerules.la \
	genericgeometry/libgenericgeometry.la \
	$(DUNE_LIBS)

include $(top_srcdir)/am/global-rules
