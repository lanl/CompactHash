########### global settings ###############
set(H_SRCS Bounds1d.h Bounds2d.h Globals1d.h Globals2d.h KDTree1d.h KDTree2d.h)

set(C_SRCS Bounds1d.c Bounds2d.c KDTree1d.c KDTree2d.c)

set(kdtree_LIB_SRCS ${C_SRCS} ${H_SRCS})

########### kdtree target ###############

add_library(kdtree SHARED ${kdtree_LIB_SRCS})

set_target_properties(kdtree PROPERTIES VERSION 1.0.0 SOVERSION 2)
install(TARGETS kdtree DESTINATION lib)

########### install files ###############

#install(FILES  KDTree.h DESTINATION include)

#========== original Makefile.am contents follow ===========

#original Makefile.am contents follow:

#default: libkdtree.la
#all: libkdtree.la
#
#AM_MAKEFLAGS = -j 4
#
#DEFAULT_INCLUDES=-I. -I..
#
#lib_LTLIBRARIES = libkdtree.la
#include_HEADERS = KDTree.h
#
#libkdtree_la_SOURCES = ${C_SRCS} ${H_SRCS}
#
