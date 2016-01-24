#!/bin/sh

#./neigh2d -t 1.0 -o br kd hg hg2 hg3 holg3 hc hc2 hc3 holc3
set -v
../neigh2d -L 7 -t 0.1 -o kd hc hc2 hc3 holc3 hg hg2 hg3 holg3 > memopt.out
fgrep Size memopt.out > memopt.dat
fgrep 256 memopt.out >> memopt.dat
