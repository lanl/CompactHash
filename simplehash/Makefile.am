default: simplehash.c libsimplehash.la
all: simplehash.c libsimplehash.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libsimplehash.la
include_HEADERS = simplehash.h

H_SRCS = simplehash.h

C_SRCS = simplehash.c

libsimplehash_la_SOURCES = ${C_SRCS} ${H_SRCS}

simplehash.c: simplehashlib_kern.inc simplehash_kern.inc

simplehashlib_kern.inc:
	./embed_source.pl simplehashlib_kern.cl > simplehashlib_kern.inc

simplehash_kern.inc:
	./embed_source.pl simplehash_kern.cl > simplehash_kern.inc

EXTRA_DIST = simplehash_kern.cl simplehashlib_kern.cl embed_source.pl

clean-local:
	rm -f *.inc
