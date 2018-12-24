#!/bin/bash

set -e

mkdir -p local
mkdir -p build
cd build

mkdir -p zlib
cd zlib
CFLAGS=-fPIC ../../zlib/configure --prefix=../../local
make
make install
cd ..

rm -rf libdivsufsort
git clone https://github.com/y-256/libdivsufsort libdivsufsort
cd libdivsufsort
cmake -DCMAKE_INSTALL_PREFIX:PATH=./ -DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON ./
make
make install
cp lib/libdivsufsort64.a ../../local/lib/
cp include/divsufsort64.h ../../local/include/
cd ..

cd ..
make
