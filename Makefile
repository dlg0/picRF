F90 = gfortran
F90FLAGS = -Jmod/ -g -fbacktrace

CUDA = nvcc
CUDA_LIB = -L/home/dg6/code/cuda/lib/ -lcudart
CUDA_FLAGS = --compiler-bindir /home/dg6/code/gcc43/usr/bin -I /usr/include/c++/4.4.1/x86_64-redhat-linux/ -I /usr/include/c++/4.4.1/ -g -arch sm_13 #-deviceemu 

EXEC = xpicRF

BLAS = ${HOME}/code/goto_blas/GotoBLAS/libgoto.a -pthread
LAPACK = ${HOME}/code/lapack/lapack-3.1.1/lapack_LINUX.a 
NETCDF = -I ${HOME}/code/netcdf/netcdf_gnu64/include -L ${HOME}/code/netcdf/netcdf_gnu64/lib -lnetcdf

LIBS = ${LAPACK} ${BLAS} ${NETCDF}
INCS = -I ${HOME}/code/netcdf/netcdf_gnu64/include

BOUNDS = -fbounds-check 
WARN = -Wall

OBJS := $(patsubst src/%.f90,obj/%.o,$(wildcard src/*.f90))
CUDA_OBJS := $(patsubst src/%.cu,obj/%.o,$(wildcard src/*.cu))

.PHONY: clean

${EXEC}: ${OBJS} ${CUDA_OBJS} 
	${F90} ${F90FLAGS} -o $@ $^ ${LIBS} ${CUDA_LIB}

obj/%.o: src/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${WARN} ${INCS}

obj/%.o: src/%.cu
	${CUDA} -c ${CUDA_FLAGS} $< -o $@

clean:
	rm obj/* mod/* x*

#	module dependencies

obj/picRF.o: obj/luxury.o obj/random.o obj/timer_class.o

#obj/picRF.o: obj/pic_kernel.o


