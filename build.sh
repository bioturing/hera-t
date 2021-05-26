#!/bin/bash

set -e

mkdir -p local
mkdir -p build
cd build

mkdir -p zlib
cd zlib
LOCAL_DIR="../../local_unix"
CFLAGS=-fPIC ../../zlib/configure --prefix=${LOCAL_DIR}
make
make install
cd ..

rm -rf libdivsufsort
git clone https://github.com/y-256/libdivsufsort libdivsufsort
cd libdivsufsort
cmake -DCMAKE_INSTALL_PREFIX:PATH=./ -DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON ./
make
make install
cp lib/libdivsufsort64.a ${LOCAL_DIR}/lib/
cp include/divsufsort64.h ${LOCAL_DIR}/include/
cd ..

cd ..
make
