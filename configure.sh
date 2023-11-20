CC=/opt/homebrew/bin/gcc-13 CXX=/opt/homebrew/bin/g++-13 \
cmake \
 -DCMAKE_INSTALL_PREFIX=$HOME/opt \
 -DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib \
 -DCMAKE_BUILD_TYPE=Release \
 -Domp=ON \
 -Dhdf5=ON \
 -Dmpi=ON \
 -DBUILD_DOC=ON \
 -DCMAKE_INCLUDE_PATH="/opt/homebrew/include;$HOME/opt/include" \
 -DCMAKE_LIBRARY_PATH="/opt/homebrew/lib;$HOME/opt/lib" \
 -DCMAKE_CXX_FLAGS="-std=c++11 -O3 -Wl,-ld_classic" \
 ..



