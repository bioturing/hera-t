#!/bin/sh

cd libdivsufsort
mkdir -p build
cd build
cmake -G "Visual Studio 14 2015 Win64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../../local -DBUILD_DIVSUVSORT64=ON ../
#cmake --build . --target install --config Release
cd ../..

cd zlib
mkdir -p build 
cd  build
cmake -G "Visual Studio 14 2015 Win64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../../local ../
#cmake --build . --target install --config Release
cd ../..

msbuild.exe Hera0.1.1/Hera0.1.1.sln "-p:Configuration=Release;Platform=x64"

