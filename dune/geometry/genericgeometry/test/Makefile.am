TOPOLOGY = CubeTopology
DIMENSION = 3

check_PROGRAMS = testsubtopologies testbasicgeometry testmatrix testsimplexify
TESTS = testsubtopologies testbasicgeometry testsimplexify


AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(ALL_PKG_LDFLAGS) $(DUNEMPILIBS)
LDADD = $(ALL_PKG_LIBS) @LDADD@ $(DUNEMPILDFLAGS)

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(AM_CPPFLAGS) -DTOPOLOGY=$(TOPOLOGY) -DDIMENSION=$(DIMENSION)

testbasicgeometry_SOURCES = testbasicgeometry.cc

testmatrix_SOURCES = testmatrix.cc
testmatrix_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
testmatrix_LDADD = $(LDADD) $(GMP_LIBS)

testsimplexify_SOURCES = testsimplexify.cc

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
