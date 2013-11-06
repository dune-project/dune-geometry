# $Id: Makefile.am 2248 2006-02-10 15:30:52Z joe $

quadraturerulesdir = $(includedir)/dune/geometry/quadraturerules
quadraturerules_HEADERS =		\
	compositequadraturerule.hh 	\
	pointquadrature.hh			\
	simplexquadrature.hh		\
	tensorproductquadrature.hh		\
	gauss_imp.hh				\
	gausslobatto_imp.hh		\
	jacobi_1_0_imp.hh			\
	jacobi_2_0_imp.hh

headercheck_IGNORE = 			\
	pointquadrature.hh			\
	simplexquadrature.hh		\
	tensorproductquadrature.hh		\
	gauss_imp.hh				\
        gausslobatto_imp.hh		\
	jacobi_1_0_imp.hh			\
	jacobi_2_0_imp.hh


noinst_LTLIBRARIES = libquadraturerules.la
libquadraturerules_la_SOURCES = quadraturerules.cc gauss.cc gausslobatto.cc jacobi_1_0.cc jacobi_2_0.cc
libquadraturerules_la_LIBADD =	\
	$(DUNE_LIBS)

extra_DIST=jacobian.mac

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
