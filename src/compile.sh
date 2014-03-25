g++ -lrt -g -L/usr/lib/x86_64-linux-gnu -lz -lbz2 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -I./ -I/usr/include/ recoverFragment.cpp -o ../bin/recoverFragment

