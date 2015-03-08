g++ -lrt -g -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 recoverFragment.cpp -o ../bin/recoverFragment -I./ -I/usr/include/ -L/usr/lib/x86_64-linux-gnu -lz -lbz2 -lpthread -lboost_system -lboost_thread
