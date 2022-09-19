#!/bin/sh
git clone https://gitlab.com/libeigen/eigen.git External/eigen
git clone https://github.com/nothings/stb.git External/stb

mkdir Build
cd Build
cmake ..
make