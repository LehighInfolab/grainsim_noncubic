module load gcc/12.4.0

wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
tar xzf hdf5-1.14.5.tar.gz
cd hdf5-1.14.5

./configure \
  --prefix=$HOME/local \
  --disable-parallel \
  --enable-shared \
  --disable-fortran \
  --disable-cxx
make -j4
make install
echo $?

# both should be no
$HOME/local/bin/h5cc -showconfig | grep -i parallel 
readelf -d $HOME/local/lib/libhdf5.so | grep -i mpi   
ls $HOME/local/lib*/libhdf5.so*