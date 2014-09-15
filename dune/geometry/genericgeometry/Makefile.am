genericgeometrydir = $(includedir)/dune/geometry/genericgeometry
genericgeometry_HEADERS = \
	codimtable.hh \
	conversion.hh \
	geometrytraits.hh \
	matrixhelper.hh \
	referencedomain.hh \
	referenceelements.hh \
	subtopologies.hh \
	topologytypes.hh

SUBDIRS = test

noinst_LTLIBRARIES = libgenericgeometry.la
libgenericgeometry_la_SOURCES = referencedomain.cc subtopologies.cc
libgenericgeometry_la_LIBADD = $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
