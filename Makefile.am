bin_PROGRAMS = neigh2d # neigh2d_memopts

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = kdtree HashFactory simplehash timer

neigh2d_C_SRCS = neigh2d.c gpu.c
neigh2d_H_SRCS = gpu.h
neigh2d_SOURCES = ${neigh2d_C_SRCS} ${neigh2d_H_SRCS}
neigh2d_LDADD = kdtree/libkdtree.la HashFactory/libHashFactory.la simplehash/libsimplehash.la timer/libtimer.la

LDADD = -lm

EXTRA_DIST = neigh2d_kern.cl APACHE_LICENSE-2_0 AUTHORS README LICENSE \
  Papers.bib NOTICE README.md ReleaseAuthorization.pdf Todo libtool

kdtree/libkdtree.la:
	cd kdtree
	make

timer/libtimer.la:
	cd timer
	make

HashFactory/libHashFactory.la:
	cd HashFactory
	make

simplehash/libsimplehash.la:
	cd simplehash
	make

clean-local:
	rm -f *.inc
	rm -f memory.out

distclean-local:
	rm -f hash-v*.tar.gz
	rm -rf build
	rm -rf autom4te.cache
	rm -f hash/stamp-h? stamp-h?
