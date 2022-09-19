@echo off
git clone https://gitlab.com/libeigen/eigen.git External
git clone https://github.com/nothings/stb.git External

mkdir Build
cd Build
cmake ..
make