default: hashit.c hasholdlib.c libhashold.la
all: hashit.c hasholdlib.c libhashold.la

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libhashold.la
include_HEADERS = hashit.h hasholdlib.h

H_SRCS = hashit.h hasholdlib.h

C_SRCS = hashit.c hasholdlib.c

libhashold_la_SOURCES = ${C_SRCS} ${H_SRCS}

hasholdlib.c: hasholdlib_kern.inc hashold_kern.inc

hasholdlib_kern.inc:
	./embed_source.pl hasholdlib_kern.cl > hasholdlib_kern.inc

hashold_kern.inc:
	./embed_source.pl hashold_kern.cl > hashold_kern.inc

EXTRA_DIST = hashold_kern.cl hasholdlib_kern.cl embed_source.pl

clean-local:
	rm -f *.inc
