F90 = gfortran
F90FLAGS = -Jmod/ -g

EXEC = xpicRF

BLAS = ${HOME}/code/goto_blas/GotoBLAS/libgoto.a -pthread
LAPACK = ${HOME}/code/lapack/lapack-3.1.1/lapack_LINUX.a 
NETCDF = -I ${HOME}/code/netcdf/netcdf_gnu64/include -L ${HOME}/code/netcdf/netcdf_gnu64/lib -lnetcdf

LIBS = ${LAPACK} ${BLAS} ${NETCDF}
INCS = -I ${HOME}/code/netcdf/netcdf_gnu64/include

BOUNDS = -fbounds-check 
WARN = -Wall

OBJS := $(patsubst src/%.f90,obj/%.o,$(wildcard src/*.f90))

.PHONY: clean

${EXEC}: ${OBJS} 
	${F90} ${F90FLAGS} -o $@ $^ ${LIBS} 

obj/%.o: src/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${BOUNDS} ${WARN} ${INCS}

clean:
	rm obj/* mod/* x*

#	module dependencies

obj/picRF.o: obj/luxury.o obj/random.o obj/timer_class.o


