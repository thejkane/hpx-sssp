export HPX_INSTALL_DIR=/N/dc2/scratch/thejkane/hpx/install
make clean
rm CMakeCache.txt
rm -rf CMakeFiles/
cmake .. \
	-DHPX_DIR=$HPX_INSTALL_DIR/lib/cmake/hpx \
	-DCMAKE_BUILD_TYPE=Debug

make VERBOSE=1
