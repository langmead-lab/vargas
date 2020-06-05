#!/bin/bash

cd htslib
autoconf && autoheader && ./configure && make -j4
cd ..

mkdir build-sse && cd build-sse
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .. && make -j4
cd ..


mkdir build-avx2 && cd build-avx2
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_AVX2_GCC=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .. && make -j4
cd ..


mkdir build-avx512bw && cd build-avx512bw
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_AVX512BW_GCC=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .. && make -j4
cd ..
