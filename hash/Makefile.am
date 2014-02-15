default: hash.c CLHash_Utilities.c hash.cl libhash.la
all: hash.c CLHash_Utilities.c hash.cl libhash.la
test: hashTest

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libhash.la
include_HEADERS = hash.h
bin_PROGRAMS = hashTest

H_SRCS = hash.h CLHash_Utilities.h

C_SRCS = hash.c CLHash_Utilities.c

libhash_la_SOURCES = ${C_SRCS} ${H_SRCS}
hashTest_SOURCES = testHelper.c hashTest.c
hashTest_LDADD = libhash.la

hash.h: hash.hp hash.hm perfectHash.hm compactHash.hm
	cpp -P -CC -x c -imacros hash.hm hash.hp | indent -linux -brf -o hash.h

hash.c: hash.cp hash.cm hash.cl perfectHash.cm compactHash.cm
	cpp -P -CC -x c -imacros hash.cm hash.cp | indent -linux -brf -o hash.c
	./embed_source.pl hash.cl >> hash.c

hash.cl:	hash.clp hash.cm perfectHash.cm compactHash.cm
	cpp -P -CC -x c -imacros hash.cm hash.clp | indent -linux -brf -o hash.cl
