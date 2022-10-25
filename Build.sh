#!/bin/sh

if [ ! -d "External/eigen" ]
then
	git clone https://gitlab.com/libeigen/eigen.git External/eigen
fi

if [ ! -d "External/stb" ]
then
	git clone https://github.com/nothings/stb.git External/stb
fi

if [ ! -d "Build" ]
then
    mkdir Build
fi

cd Build
cmake ..
make
