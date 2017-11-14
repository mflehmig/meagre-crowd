# TODO

* bebop and sparse_matrix_converter built with debug symbols!


# Configure and Install

## Get Sources
    git clone https://github.com/mflehmig/meagre-crowd.git meagre-crowd.git

## Build Dependencies

The script get_dependencies.sh automaticaly downloads and builds the dependencies. It did not work for me. I would have make do much adapations, thus, I built the dependencies by my own.

    cd meagre-crowd.git
    mkdir 3rdParty && cd 3rdParty
    wget http://bebop.cs.berkeley.edu/smc/tarballs/bebop_make.tar.gz
    wget http://bebop.cs.berkeley.edu/smc/tarballs/bebop_util.tar.gz
    wget http://bebop.cs.berkeley.edu/smc/tarballs/sparse_matrix_converter.tar.gz

    tar xvfz bebop_make.tar.gz
    tar xvfz bebop_util.tar.gz
    tar xvfz bebop_sparse_matrix_converter.tar.gz

Adopt *bebop_make/Makefile.include.linux* to your needs/wishes/requirements, e.g., *CC = icc*.

    cd bebop_util/ && make

    cd ..
    mkdir lib
    mkdir -p include/bebop/smc
    cd sparse_matrix_converter/ && make && cd ../

    cp sparse_matrix_converter/include/bebop/smc/*.h include/bebop/scm/
    cp sparse_matrix_converter/lib* lib/
    cp bebop_util/include/bebop/util/enumerations.h include/bebop/util/
    cp bebop_util/include/bebop/util/init.h include/bebop/util/
    cp bebop_util/libbebop_util.* lib/


## Configure and Build meagre-crowd

    cd ..
    ./autogen
    ./configure CC=mpicc CFLAGS='/PATH/meagre-crowd.git/3rdParty/include/' \
                LDFLAGS='-L/PATH/meagre-crowd.git/3rdParty/lib/ -Wl,-rpath,/PATH/meagre-crowd.git/3rdParty/lib/' \
                CPPFLAGS="-I/PATH/meagre-crowd.git/3rdParty/include/ -I/usr/include/suitesparse/ -I/usr/include/openmpi-x86_64/superlu_dist/ \
                [--prefix=/PATH/TO/INSTALLATION]
    make
    [make install]


### Taurus
  module load intelmpi/2017.2.174
  module load mkl/2017
  
  module load gcc/7.1.0 openmpi/2.1.0-gnu6.3


  mkdir 3rdParty && cd 3rdParty
  wget http://bebop.cs.berkeley.edu/smc/tarballs/bebop_make.tar.gz
  wget http://bebop.cs.berkeley.edu/smc/tarballs/bebop_util.tar.gz
  wget http://bebop.cs.berkeley.edu/smc/tarballs/sparse_matrix_converter.tar.gz

  tar xvfz bebop_make.tar.gz
  tar xvfz bebop_util.tar.gz
  tar xvfz bebop_sparse_matrix_converter.tar.gz
  In bebop_make/Makefile.include.linux:
    CC = icc
    LINKER = icc
    
  cd bebop_util/ && make
  
  mkdir lib
  mkdir -p include/bebop/scm
  
  cp sparse_matrix_converter/include/bebop/smc/*.h include/bebop/scm/
  cp sparse_matrix_converter/lib* lib/
  cp bebop_util/include/bebop/util/enumerations.h include/bebop/util/
  cp bebop_util/include/bebop/util/init.h include/bebop/util/
  cp bebop_util/libbebop_util.* lib/


## Configure and Build meagre-crowd
  cd ..
    export LIBRARY_PATH=$LIBRARY_PATH:/home/mflehmig/Projekte/PARADOM/meagre-crowd.git/3rdParty/lib/
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mflehmig/Projekte/PARADOM/meagre-crowd.git/3rdParty/lib/ (reicht nur LIBRARY_PATH?)
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/home/mflehmig/Projekte/PARADOM/meagre-crowd.git/3rdParty/include/
    
    // Enable SuperLU_DIST
    export LIBRARY_PATH=$LIBRARY_PATH:$SUPERLU_DIST_LIB
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:$SUPERLU_DIST_INC

    // Enable UMFPACK and cholmod
    export LIBRARY_PATH=$LIBRARY_PATH:$SUITESPARSE_LIB
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:$SUITESPARSE_INC


    ./autogen
    ./configure CC=mpicc
    make









# Known Issues and Errors

## Runtime
* input error: Failed to load matrix:
--> Solution: Adopt matrix file, like
  1 1 .72827 --> 1 1 0.72827
  1 1 -.72827 --> 1 1 -0.72827


## Configuration and Installation

* bebop.lib not found:
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mf/opt/PARADOM/AP3/meagre-crowd-3rdParty/bebop_util/:/home/mf/opt/PARADOM/AP3/meagre-crowd-3rdParty/sparse_matrix_converter/

* mpi.h not found:
  export C_INCLUDE_PATH=$C_INCLUDE_PATH//usr/include/openmpi-x86_64/


# Initial Adaptions to Build the Code

* inline Keyword bei Funktion (siehe https://gcc.gnu.org/gcc-5/porting_to.html)
  file.h:
    /*inline*/ unsigned int matrix_rows( const matrix_t* const A );

  solvers.h
    /*inline*/ int solver_uses_mpi( const int solver );
    /*inline*/ int solver_requires_mpi( const int solver );
    /*inline*/ int solver_uses_omp( const int solver );
    /*inline*/ int solver_requires_omp( const int solver );

* Use of Matio datatypes:
  file.c
    sparse_t* st = t->data; --> mat_sparse_t* st = t->data;



src/matrix.o: In function `copy_matrix':
matrix.c:(.text+0x536): undefined reference to `_data_width'
src/matrix.o: In function `cmp_matrix':
matrix.c:(.text+0xc66): undefined reference to `_data_width'
src/matrix.o: In function `_coo2drow':
matrix.c:(.text+0x1a42): undefined reference to `_data_width'
src/matrix.o: In function `_drow2coo':
matrix.c:(.text+0x1bf7): undefined reference to `_data_width'
src/matrix.o: In function `_drow2dcol':
matrix.c:(.text+0x1eaf): undefined reference to `_data_width'
src/matrix.o:matrix.c:(.text+0x1eea): more undefined references to `_data_width' follow
collect2: error: ld returned 1 exit status


solver_superlu_dist.c:
  typedef struct {
    superlu_options_t options;  --> superlu_dist_options_t options;

  void solver_evaluate_superlu_dist( solver_state_t* s, matrix_t* b, matrix_t* x ) {
    ...
    LUstructInit(p->A.nrow, p->A.ncol, &(p->lu)); --> LUstructInit(p->A.ncol, &(p->lu));

Moved size_t _data_width( const enum matrix_data_type_t t ); from matrix.c to matrix.h





