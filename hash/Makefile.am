default: hash.c CLHash_Utilities.c hash.cl libhash.la
all: hash.c CLHash_Utilities.c hash.cl libhash.la
test: hashTest

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libhash.la
include_HEADERS = hash.h
bin_PROGRAMS = hashTest

hashTest_CPPFLAGS = -DHAVE_OPENCL -DHAVE_CL_DOUBLE -I$(OPENCL_INCLUDE)

H_SRCS = hash.h CLHash_Utilities.h

C_SRCS = hash.c CLHash_Utilities.c

libhash_la_SOURCES = ${C_SRCS} ${H_SRCS}
hashTest_SOURCES = testHelper.c hashTest.c
hashTest_LDADD = libhash.la

hash.h: hash.hp hash.hm perfectHash.hm compactHash.hm
	${CPREPROCESSOR} -CC -x c -imacros hash.hm hash.hp | ${INDENT} -linux -brf -o hash.h

hash.c: hash.cp hash.cm hash.cl perfectHash.cm compactHash.cm
	${CPREPROCESSOR} -CC -x c -imacros hash.cm hash.cp | ${INDENT} -linux -brf -o hash.c
	./embed_source.pl hash.cl >> hash.c

hash.cl:	hash.clp hash.cm perfectHash.cm compactHash.cm
	${CPREPROCESSOR} -CC -x c -imacros hash.cm hash.clp | ${INDENT} -linux -brf -o hash.cl

EXTRA_DIST = perfectHash.cm perfectHash.hm primetest.c testHelper.h hash.hp \
    hash.hm hash.cp hash.cm hash.h hash.c hash.cl hash.clp embed_source.pl \
    compactHash.cm compactHash.hm baseHash.cm baseHash.hm
