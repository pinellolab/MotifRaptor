#!/bin/bash
cd MotifRaptor/thirdparty/libdivsufsort
cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="./"
make
make install
cd examples
make
cp mksary ../../../SNPScanner

