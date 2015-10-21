default: simplehash.c simplehashlib.c libsimplehash.la
all: simplehash.c simplehashlib.c libsimplehash.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libsimplehash.la
include_HEADERS = simplehash.h simplehashlib.h

H_SRCS = simplehash.h simplehashlib.h

C_SRCS = simplehash.c simplehashlib.c

libsimplehash_la_SOURCES = ${C_SRCS} ${H_SRCS}

simplehashlib.c: simplehashlib_kern.inc simplehash_kern.inc

simplehashlib_kern.inc:
	./embed_source.pl simplehashlib_kern.cl > simplehashlib_kern.inc

simplehash_kern.inc:
	./embed_source.pl simplehash_kern.cl > simplehash_kern.inc

EXTRA_DIST = simplehash_kern.cl simplehashlib_kern.cl embed_source.pl

clean-local:
	rm -f *.inc
