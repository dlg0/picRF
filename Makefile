F90 = gfortran

EXEC = xpicRF
BLAS = /home/dg6/code/goto_blas/GotoBLAS/libgoto.a -pthread
LAPACK = /home/dg6/code/lapack/lapack-3.1.1/lapack_LINUX.a 
NETCDF = -I /home/dg6/code/netcdf/netcdf_gnu64/include -L /home/dg6/code/netcdf/netcdf_gnu64/lib -lnetcdf

OBJS = 
BOUNDS = -fbounds-check -g
WARN = -Wall

${EXEC}: ${OBJS} src/picRF.f90
	${F90} src/picRF.f90 -o ${EXEC} ${BOUNDS} ${LAPACK} ${BLAS} ${WARN} ${NETCDF}

clean:
	rm obj/* mod/* x*
