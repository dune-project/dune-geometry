genericgeometrydir = $(includedir)/dune/geometry/genericgeometry
genericgeometry_HEADERS = cachedmapping.hh  codimtable.hh  conversion.hh \
                          cornermapping.hh  geometrytraits.hh \
                          geometry.hh  hybridmapping.hh \
                          mapping.hh  mappingprovider.hh  matrixhelper.hh \
                          maximum.hh  referencedomain.hh  referenceelements.hh \
                          subtopologies.hh \
                          topologytypes.hh  traceprovider.hh

SUBDIRS = test

noinst_LTLIBRARIES = libgenericgeometry.la
libgenericgeometry_la_SOURCES = referencedomain.cc subtopologies.cc
libgenericgeometry_la_LIBADD = $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
