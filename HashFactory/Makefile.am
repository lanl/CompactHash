default: HashFactory.c CLHash_Utilities.c HashFactory.cl libHashFactory.la
all: HashFactory.c CLHash_Utilities.c HashFactory.cl libHashFactory.la
test: hashTest

AM_MAKEFLAGS = -j 4

DEFAULT_INCLUDES=-I. -I..

lib_LTLIBRARIES = libHashFactory.la
include_HEADERS = HashFactory.h
bin_PROGRAMS = hashTest

hashTest_CPPFLAGS = -DHAVE_OPENCL -DHAVE_CL_DOUBLE -I$(OPENCL_INCLUDE)

H_SRCS = HashFactory.h CLHash_Utilities.h

C_SRCS = HashFactory.c CLHash_Utilities.c

libHashFactory_la_SOURCES = ${C_SRCS} ${H_SRCS}
hashTest_SOURCES = testHelper.c hashTest.c
hashTest_LDADD = libHashFactory.la

HashFactory.h: HashFactory.hp HashFactory.hm perfectHash.hm compactHash.hm
	${CPREPROCESSOR} -CC -x c -imacros HashFactory.hm HashFactory.hp | ${INDENT} -linux -brf -o HashFactory.h

HashFactory.c: HashFactory.cp HashFactory.cm HashFactory.cl perfectHash.cm compactHash.cm HashFactory.h
	${CPREPROCESSOR} -CC -x c -imacros HashFactory.cm HashFactory.cp | ${INDENT} -linux -brf | sed -e '/OMP_PRAGMA_PARALLEL_FOR/s/[[:space:]]*OMP_PRAGMA_PARALLEL_FOR;/#pragma omp parallel for/' > HashFactory.c
	./embed_source.pl HashFactory.cl >> HashFactory.c

HashFactory.cl:	HashFactory.clp HashFactory.cm perfectHash.cm compactHash.cm
	${CPREPROCESSOR} -CC -x c -imacros HashFactory.cm HashFactory.clp | ${INDENT} -linux -brf -o HashFactory.cl

EXTRA_DIST = perfectHash.cm perfectHash.hm primetest.c testHelper.h HashFactory.hp \
    HashFactory.hm HashFactory.cp HashFactory.cm HashFactory.h HashFactory.c HashFactory.cl HashFactory.clp embed_source.pl \
    compactHash.cm compactHash.hm baseHash.cm baseHash.hm
