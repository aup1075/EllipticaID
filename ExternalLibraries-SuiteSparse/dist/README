#===========================================================
# H O W   TO   C O M P I L E   U M F P A C K   L I B R A R Y
#===========================================================

# General note:
# intel compilers:
#       for C       = icc
#       optimization flag for icc = -Ofast
#       for fortran = ifort
#       optimiza flag for ifort   = -fast

# follow the commands below step by step:

####################
## 1. Compile Lapack
####################
$ tar -zxvf v3.9.0.tar.gz
$ cd lapack-3.9.0
$ cp make.inc.example make.inc
$ edit make.inc
  # a. pick your compiler	
  # b. edit FFLAGS: -O2 -> -O3
  # c. add `-fPIC` to CFLAGS and FFLAGS
  # d. edit BLASLIB: $(TOPSRCDIR)/librefblas.a -> $(TOPSRCDIR)/libblas.a
$ make -j

# Notes:
# there might be some errors about timer during compile time, 
# which one can adjust the TIMER flag in "make.inc". For instance,
# set TIMER = INT_CPU_TIME when using ifort compiler, or TIMER = NONE.

# Note:
# one can also adjust Makefile, for example, removing blas_testing or 
# lapack_testing.

#####################
## 2. Compile UMFPACK
#####################
$ tar -zxvf SuiteSparse-5.7.2.tar.gz
$ cd SuiteSparse-5.7.2
$ cd SuiteSparse_config
$ edit SuiteSparse_config.mk
  # go to 'optimization' and adjust your optimization flag
  # go to `required libraries` section:
  # a. add -lgfortran to `LDLIBS` for fortran compiler;
  #    for intel compiler add "-ldl -lirc -lifcore -ldl" to `LDLIBS`. 
  # b. add BLAS = -L/path/to/lapack -lblas
  # c. add LAPACK=-L/path/to/lapack -llapack
  # go to UMFPACK configuration section:
  # a. remove the comment for `UMFPACK_CONFIG = -DNCHOLMOD`
  #    this makes the UMFPACK be dependent only on AMD

$ make
$ cd ../ && cd AMD && make
$ cd ../ && cd UMFPACK && make

# you can find the `libumfpack.a` in uiteSparse-5.7.2/Lib to link.


#=====================================================
# H O W   TO   C O M P I L E   S I L O   L I B R A R Y
#=====================================================

# follow the commands below step by step:

###################
## 1. Compile zlib:
###################
$ tar -zxvf zlib-1.2.11.tar.gz
$ cd zlib-1.2.11
$ ./configure --prefix=/absolute/path/to/zlib/dir
$ make test
$ make install

###################
## 2. Compile hdf5:
###################
$ tar -zxvf hdf5-1.12.0.tar.gz
$ cd hdf5-1.12.0
$ mkdir HDF5_compiled # this is required otherwise bin dir conflicts
$ ./configure --prefix=/absolute/path/to/HDF5_compiled \
  	      --with-zlib=path/to/zlib/directory       \
              --enable-build-mode=production           \
              --with-default-api-version=v110          \
              --enable-threadsafe --disable-hl         \
              --enable-optimization=high               \
              --disable-internal-debug

$ ./configure --help # for more help if you have spacial case
$ make -j               # paralled compilation is recommended
$ make -j check         # to doc test it's not mandatory, skip it
$ make -j install       # to install
$ make -j check-install # to doc test the installation

###################
## 3. Compile silo:
###################
$ tar -zxvf silo-4.10.2.tar.gz
$ cd silo-4.10.2
$ mkdir silo_compiled
$ ./configure --prefix=/absolute/path/to/silo_compiled \
              CFLAGS=-O3 --enable-optimization         \
   --with-hdf5=/path/to/HDF5_compiled/include,/path/to/HDF5_compiled/lib

$ ./configure --help # for more help if you have spacial case

$ make -j
$ make -j install




       



