# $Id: Makefile.am 6217 2010-01-20 23:02:03Z christi $

if BUILD_DOCS
PAGES =
endif

SUBDIRS = refelements

# place in website
BASEDIR=../..
CURDIR=doc/appl

# install the html pages
docdir=$(datadir)/doc/dune-geometry/appl
if BUILD_DOCS
DOCFILES = $(PAGES)
EXTRA_DIST = CMakeLists.txt $(PAGES)
endif

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
