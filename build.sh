export HPX_INSTALL_DIR=/Users/thejaka/classes/quals/andrew/hpx/install
make clean
#rm CMakeCache.txt
cmake -DHPX_DIR=$HPX_INSTALL_DIR/lib/cmake/hpx .
make VERBOSE=1
