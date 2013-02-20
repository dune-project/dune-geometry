# $Id: duneproject 6399 2011-03-18 16:05:35Z mblatt $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST= CMakeLists.txt dune.module

SUBDIRS = m4 dune doc lib

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

