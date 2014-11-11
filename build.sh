#export HPX_INSTALL_DIR=/Users/thejaka/classes/quals/andrew/hpx/install
#export HPX_INSTALL_DIR=/Users/thejaka/classes/quals/andrew/hpx-latest/install
export HPX_INSTALL_DIR=/Users/thejaka/classes/quals/andrew/hpx/hpx/build
make clean
rm CMakeCache.txt
rm -rf CMakeFiles
cmake -DHPX_DIR=$HPX_INSTALL_DIR/lib/cmake/hpx -DCMAKE_BUILD_TYPE=Debug ..
make VERBOSE=1
