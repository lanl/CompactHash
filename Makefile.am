bin_PROGRAMS = neigh2d # neigh2d_memopts

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = kdtree hash oldHash timer

neigh2d_C_SRCS = neigh2d.c gpu.c
neigh2d_H_SRCS = gpu.h
neigh2d_SOURCES = ${neigh2d_C_SRCS} ${neigh2d_H_SRCS}
neigh2d_LDADD = kdtree/libkdtree.la hash/libhash.la oldHash/libhashold.la timer/libtimer.la

LDADD = -lm

EXTRA_DIST = neigh2d_kern.cl config.h.in APACHE_LICENSE-2_0 AUTHORS README LICENSE Papers.bib NOTICE

kdtree/libkdtree.la:
	cd kdtree
	make

timer/libtimer.la:
	cd timer
	make

hash/libhash.la:
	cd hash
	make

oldHash/libhashold.la:
	cd oldHash
	make

clean-local:
	rm -f *.inc
	rm -f memory.out

distclean-local:
	rm -f config.h.in~
	rm -f hash-v*.tar.gz
	rm -rf build
	rm -rf autom4te.cache
	rm -f hash/stamp-h? stamp-h?
