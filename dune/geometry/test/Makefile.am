# which tests to run
TESTS = test-axisalignedcubegeometry test-quadrature test-genericquadrature test-referenceelements test-high-dim-refelements $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CXXFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs
test_axisalignedcubegeometry_SOURCES = test-axisalignedcubegeometry.cc

test_referenceelements_SOURCES = test-referenceelements.cc

test_high_dim_refelements_SOURCES = test-high-dim-refelements.cc

test_quadrature_SOURCES = test-quadrature.cc

test_genericquadrature_SOURCES = test-genericquadrature.cc

include $(top_srcdir)/am/global-rules
