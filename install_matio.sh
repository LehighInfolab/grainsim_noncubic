module load gcc/12.4.0 cmake/3.30.5 # hdf5
cmake --version

wget https://github.com/tbeu/matio/releases/download/v1.5.28/matio-1.5.28.tar.gz
tar xzf matio-1.5.28.tar.gz
cd matio-1.5.28
rm -rf build && mkdir build && cd build

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/local \
  -DMATIO_MAT73=ON \
  -DMATIO_SHARED=ON \
  -DBUILD_TESTING=OFF \
  -DHDF5_ROOT=$HOME/local \
  -DHDF5_PREFER_PARALLEL=OFF
make -j4
make install
echo $?

ls $HOME/local/include/matio.h
ls $HOME/local/lib*/libmatio.so*