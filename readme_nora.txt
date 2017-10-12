./configure CC=mpicc CPPFLAGS="-I/home/hundertsieben//Documents/meagre-crowd/3rdParty/include/ -I/usr/lib/suitesparse/include/ -I/home/hundertsieben/opt/SuperLU_DIST_5.2.0/include/" LDFLAGS="-L/usr/lib/suitesparse/lib/ -Wl,-rpath,/usr/lib/suitesparse/lib/ -L/home/hundertsieben/opt/SuperLU_DIST_5.2.0/lib/ -L/home/hundertsieben/Documents/meagre-crowd/3rdParty/lib/ -Wl,-rpath,/home/hundertsieben/Documents/meagre-crowd/3rdParty/lib/ -lbebop_util -L/usr/lib/x86_64-linux-gnu/ -lmetis -L/usr/lib/ -lblas"

export LD_LIBRARY_PATH=./3rdParty/lib/
