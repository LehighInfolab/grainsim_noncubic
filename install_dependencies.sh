#!/bin/bash
set -euo pipefail

PREFIX=$HOME/local
JOBS=$(nproc)
module load gcc/12.4.0 cmake/3.30.5

fetch() {  # fetch <name> <version> <url>
    wget -O "$1-$2.tar.gz" "$3"
    rm -rf "$1-$2"
    tar xzf "$1-$2.tar.gz"
}

build_libzip() {
    [ -f "$PREFIX/lib/libzip.so" ] && { echo "libzip: already installed, skipping"; return; }
    fetch libzip 1.10.1 https://libzip.org/download/libzip-1.10.1.tar.gz
    cd libzip-1.10.1 && rm -rf build && mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_INSTALL_LIBDIR=lib \
        -DENABLE_CRYPTO=OFF -DENABLE_OPENSSL=OFF \
        -DENABLE_GNUTLS=OFF -DENABLE_MBEDTLS=OFF
    make -j$JOBS && make install
    cd ../..
}

build_hdf5() {
    [ -f "$PREFIX/lib/libhdf5.so" ] && { echo "hdf5: already installed, skipping"; return; }
    fetch hdf5 1.14.5 https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
    cd hdf5-1.14.5
    ./configure --prefix=$PREFIX --disable-parallel --enable-shared \
        --disable-fortran --disable-cxx
    make -j$JOBS && make install
    cd ..
}

build_matio() {
    [ -f "$PREFIX/lib/libmatio.so" ] && { echo "matio: already installed, skipping"; return; }
    fetch matio 1.5.28 https://github.com/tbeu/matio/releases/download/v1.5.28/matio-1.5.28.tar.gz
    cd matio-1.5.28 && rm -rf build && mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_INSTALL_LIBDIR=lib \
        -DMATIO_MAT73=ON -DMATIO_SHARED=ON -DBUILD_TESTING=OFF \
        -DHDF5_ROOT=$PREFIX -DHDF5_PREFER_PARALLEL=OFF \
        -DCMAKE_INSTALL_RPATH=$PREFIX/lib -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON
    make -j$JOBS && make install
    cd ../..
}

verify() {
    echo "=== verify ==="
    ldd $PREFIX/lib/libzip.so   | grep -qiE 'crypto|ssl' && { echo "FAIL: libzip links openssl"; exit 1; }
    ldd $PREFIX/lib/libmatio.so | grep -qi  mpi           && { echo "FAIL: matio chain links MPI"; exit 1; }
    ldd $PREFIX/lib/libmatio.so | grep hdf5 | grep -q "$PREFIX" || { echo "FAIL: matio using wrong hdf5"; exit 1; }
    $PREFIX/bin/h5cc -showconfig | grep -i 'parallel hdf5'
    echo "Installed in: $PREFIX/lib"
}

mkdir -p "$HOME/build-deps" && cd "$HOME/build-deps"
build_libzip
build_hdf5
build_matio
verify