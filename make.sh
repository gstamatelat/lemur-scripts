#!/bin/bash

OUT=rel
LIBS="-lindri -lz -lpthread -lm"
INCLUDE="/home/james/work/scrambler/indri-5.5-src/include~/work/scrambler/indri-5.5-src/include/"
LIB=""

#g++ -O3 -Wall -DP_NEEDS_GNU_CXX_NAMESPACE=1 rel.cpp\
# -fopenmp -D_GLIBCXX_PARALLEL\
# -march=native\
# -o rel -I$INCLUDE -L$LIB $LIBS

#g++ -O3 -Wall -DP_NEEDS_GNU_CXX_NAMESPACE=1 stem.cpp\
# -fopenmp -D_GLIBCXX_PARALLEL\
# -march=native\
# -o stem -I$INCLUDE -L$LIB $LIBS

#g++ -O3 -Wall -DP_NEEDS_GNU_CXX_NAMESPACE=1 sim.cpp\
# -fopenmp -D_GLIBCXX_PARALLEL\
# -march=native\
# -o sim -I$INCLUDE -L$LIB $LIBS

#g++ -O3 -Wall -DP_NEEDS_GNU_CXX_NAMESPACE=1 rel-cover.cpp\
# -fopenmp -D_GLIBCXX_PARALLEL\
# -march=native\
# -o rel-cover -I$INCLUDE -L$LIB $LIBS

g++ -O3 -Wall -DP_NEEDS_GNU_CXX_NAMESPACE=1 tokenize.cpp\
 -fopenmp -D_GLIBCXX_PARALLEL\
 -march=native\
 -o tokenize -I$INCLUDE -L$LIB $LIBS
