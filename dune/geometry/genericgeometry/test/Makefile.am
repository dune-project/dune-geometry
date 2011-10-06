GEOMETRYTYPE = cube
DIMENSION = 3

check_PROGRAMS = testbasicgeometry testsubtopologies testmatrix
TESTS = testbasicgeometry testsubtopologies

AM_CPPFLAGS = @AM_CPPFLAGS@ $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(ALL_PKG_LDFLAGS) $(DUNEMPILIBS)
LDADD = $(ALL_PKG_LIBS) @LDADD@ $(DUNEMPILDFLAGS)

testbasicgeometry_SOURCES = testbasicgeometry.cc

testsubtopologies_SOURCES = testsubtopologies.cc
testsubtopologies_CPPFLAGS = $(AM_CPPFLAGS)			\
	-DGEOMETRYTYPE=$(GEOMETRYTYPE) -DDIMENSION=$(DIMENSION)

testmatrix_SOURCES = testmatrix.cc
testmatrix_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
testmatrix_LDADD = $(LDADD) $(GMP_LIBS)

include $(top_srcdir)/am/global-rules
