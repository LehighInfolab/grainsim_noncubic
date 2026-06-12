# load cmake 
module load gcc/12.4.0 cmake/3.30.5
cmake --version

wget https://libzip.org/download/libzip-1.10.1.tar.gz
tar xzf libzip-1.10.1.tar.gz
cd libzip-1.10.1
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local
make -j4
make install
echo $? 

ls $HOME/local/include/zip.h
ls $HOME/local/lib64/libzip.so*