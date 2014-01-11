
# $Id: Makefile.am 7400 2011-02-09 21:27:22Z christi $

# include this into website
EXTRAINSTALL=gg_hexahedron.png gg_line.png gg_prism.png gg_pyramid.png \
  gg_quadrilateral.png gg_tetrahedron.png gg_triangle.png\
  gg_hexahedron_edges.png gg_tetrahedron_edges.png gg_prism_edges.png gg_pyramid_edges.png

# place in website
BASEDIR=../../..
CURDIR=doc/appl/refelements

# install docs
docdir=$(datadir)/doc/dune-geometry/appl/refelements
if BUILD_DOCS
EXTRA_DIST= CMakeLists.txt $(EXTRAINSTALL)
DOCFILES= $(EXTRAINSTALL)
endif

if INKSCAPE
update-images: $(EXTRAINSTALL)
endif

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/inkscape.am
include $(top_srcdir)/am/global-rules
