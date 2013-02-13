
# $Id: Makefile.am 7400 2011-02-09 21:27:22Z christi $

# include this into website
EXTRAINSTALL=hexahedron.png line.png prism.png pyramid.png \
  quadrilateral.png tetrahedron.png triangle.png\
  hexahedron_edges.png tetrahedron_edges.png prism_edges.png pyramid_edges.png

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
