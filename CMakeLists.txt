# set up project
project("dune-geometry" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.8)

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("lib")

export(TARGETS dunegeometry FILE DuneGeometryTargets.cmake)
# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
