make clean
./configure --disable-hslab --disable-hlmtools --enable-hdecode --prefix=$PWD
make
make install
