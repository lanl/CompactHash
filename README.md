Compact Hash Algorithms for Computational Meshes: Neighbors 2D Testbed 
======================================================================

This project explores compact hash algorithms for CPUs and GPUs using 
OpenCL. This project contains compact hash libraries and a performance testing
infrastructure. The libraries are:

This code is a set of hash functions applied to a sample mesh. It supports the 
paper "Compact Hash Algorithms for Computational Meshes," submitted to the 
SIAM Journal of Scientific Computing.

hash/libhash.la
===============

This is a hash library with a hash factory that allows multiple hash tables
to exist at one time. The library supports perfect and compact hash methods
based on different algorithms. The code is auto generated for CPU and GPUs
in the OpenCL language for int32 and int64 datatypes from a single source.
The auto generation uses cpp macros so that it is compatible for C++, C and
Fortran applications.

oldHash/libhashold.la
=====================

A simpler hash library that autoselects between perfect and compact hash
methods. This library supports CPU and GPUs using OpenCL. Only one
hash table at a time is allowed. The algorithms are a similar, but simpler
version of those in libhash.la. The same OpenCL kernel and function
embedding is used in both libraries. It will be easier to understand the
source code and algorithms in this version of the library since it does not
have code generation using macros.

See the README for instructions on building using configure and 
running a variety of performance tests. 

Authors:
========

Bob Robey XCP-2 (brobey@lanl.gov)
Peter Ahrens XCP-2 (peter.ahrens@lanl.gov, ptrahrens@gmail.com)
Sara Hartse XCP-2 (sara@lanl.gov, sara.hartse@gmail.com)
Rebecka Tumblin (rtumblin@lanl.gov, rebeckatumblin@gmail.com)





 
