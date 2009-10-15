#include <stdio.h>
#include "cuda.h"
#include <culapackdevice.h> 


#define CUERR { cudaError_t err; \
	if (( err = cudaGetLastError()) != cudaSuccess ) \
	{ \
		printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); \
		printf("Thread aborting ...\n"); \
		return; \
	} }

#define NHIST 256
#define NWARP 6
#define NBLOCK 64 
#define NHISTNWARP (NHIST * NWARP)
#define NHISTNBLOCK (NHIST * NBLOCK)
#ifdef  __DEVICE_EMULATION__
	#define WARP_LOG2SIZE 0
#else
	#define WARP_LOG2SIZE 5
#endif
#define NTHREAD_PB (NWARP << WARP_LOG2SIZE) 

#define SMEM_ATOMICS
#define GMEM_ATOMICS

__global__ void initialize(float *u, int n)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if( i==0 ) u[i] = 0.0;
    else if( i<n ) u[i] = 0.0;
}

__global__
void compute_update(float *u, float *du, int n)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if( i>0 && i<n-1 )
    {   
    float u_left = u[i-1],
    u_i = u[i],
    u_right = u[i+1];
    du[i] = (u_left + u_right)/2 - u_i;
    }   
}


__global__
void apply_update(float *u, float lambda,float *du, int n)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if( i>0 && i<n-1 )
    u[i] = u[i] + lambda*du[i];
}


__global__ void remove_average ( float *data, float dx ){

	__shared__ float total;
	__shared__ float e0;

	e0	= 8.85418782e-12;

	total	= 0.0;
	for (int i=0;i<NHIST;i++){
		total += data[i];
	}
	for (int i=0;i<NHIST;i++){
		data[i] = -(data[i] - total / NHIST) / e0 * dx * dx;
	}

	//printf("total: %9.7f\n", total);	

}

__global__ void compute_eField ( float *phi, float *E, float dx ){

	E[0]	= ( phi[NHIST-1] - phi[1] ) / ( 2 * dx );
	E[NHIST-1]	= ( phi[NHIST-2] - phi[0] ) / ( 2 * dx ); 

	for ( int i=1;i<NHIST-2;i++ ){

		E[i]	= ( phi[i-1] - phi[i+1] ) / ( 2 * dx );

	}

}

#ifdef SMEM_ATOMICS
	
	__device__ inline void warp_hist_kernel ( unsigned int *s_warp_hist, unsigned int dataIdx, unsigned int threadTag )
	{
		atomicAdd ( s_warp_hist + dataIdx, 1 );
	}
	
#else

	__device__ inline void warp_hist_kernel ( volatile unsigned int *s_warp_hist, unsigned int dataIdx, unsigned int threadTag )
	{
	
		unsigned int count;
		do{
			count = s_warp_hist[dataIdx] & 0x07FFFFFFU;
			count = threadTag | ( count + 1 );
			s_warp_hist[dataIdx] = count;
		}while(s_warp_hist[dataIdx] != count );
		
	}

#endif

static __global__ void hist_kernel ( float *xIn, int nx, unsigned int *histOutBlock, float xMin, float xRng, float *rho, float weight, float dx)
{
	const int myThrdId = blockIdx.x * blockDim.x + threadIdx.x;
	const int nThreads = blockDim.x * gridDim.x;

	const float e = 1.602e-19;
	unsigned int ii;
	float inv_xRng;
	float inv_dx;

	inv_xRng	= 1.0 / xRng;
	inv_dx	= 1.0 / dx;

	// create shared memory copy of histogram for each warp
	// but there a 6 warps in the block hence 6*nhist
	#ifdef SMEM_ATOMICS
		__shared__ unsigned int s_warp_hist[NHISTNWARP];
	#else
		volatile __shared__ unsigned int s_warp_hist[NHISTNWARP];
	#endif

	// intitialise these 

	for (int i=threadIdx.x;i<NHISTNWARP;i+=blockDim.x)
	{
		s_warp_hist[i]	= 0;
	}


	#ifdef __DEVICE_EMULATION__
		const unsigned int threadTag = 0;
	#else
		const unsigned int threadTag = threadIdx.x << (32 - WARP_LOG2SIZE);
	#endif

	const int warpOffset = ( threadIdx.x >> WARP_LOG2SIZE ) * NHIST;

	__syncthreads();

	// work only on my chunk of the particle list
	for (int i = myThrdId; i < nx; i += nThreads)
	{
		
		// calculate histogram index for this data
		ii = rint((xIn[i]-xMin) * inv_xRng * NHIST);

		// mod (%) operation here is for periodic boundary conditions, 
		// i.e., 257 -> 1 
		warp_hist_kernel ( s_warp_hist+warpOffset, ii%256, threadTag );

	}

	__syncthreads();	

	// reduce the per warp histograms to a per block one
	for (int i=threadIdx.x; i<NHIST; i+=blockDim.x)
	{
		unsigned int warp_sum = 0;

		for (int j = 0; j < NHISTNWARP; j += NHIST)
		{
			warp_sum += s_warp_hist[j+i] & 0x07FFFFFFU;
		}

		// copy to per block global histogram

		#ifdef GMEM_ATOMICS
			atomicAdd ( histOutBlock + i, warp_sum );
		#else
			histOutBlock[blockIdx.x * NHIST + i]	= warp_sum;
		#endif
	}

	// convert histogram to real units of charge density

	for (int i=threadIdx.x; i<NHIST; i+=blockDim.x)
	{
		rho[i]	= (-1.0) * e * weight * inv_dx * histOutBlock[i];
	}
	

}

#ifndef GMEM_ATOMICS
	static __global__ void mergePerBlockHistograms_kernel ( unsigned int *d_Result ){
	
	    __shared__ unsigned int data[NBLOCK];
	
	    data[threadIdx.x] = d_Result[threadIdx.x * NHIST + blockIdx.x];
	
	    for(int stride = NBLOCK / 2; stride > 0; stride >>= 1){
	        __syncthreads();
	        if(threadIdx.x < stride)
	            data[threadIdx.x] += data[threadIdx.x + stride];
	    }
	
	    if(threadIdx.x == 0)
	        d_Result[blockIdx.x] = data[0];
	}
#endif

void checkStatus(culaStatus status)
{
    if(!status)
        return;

    if(status == culaArgumentError)
        printf("Invalid value for parameter %d\n", culaGetErrorInfo());
    else if(status == culaRuntimeError)
        printf("Runtime error (%d)\n", culaGetErrorInfo());
    else
        printf("%s\n", culaGetStatusString(status));

    culaShutdown();
    exit(EXIT_FAILURE);
}

/*
	fortran calls this subroutine with
	arguments passed by reference
*/

extern "C" void cudahist_ ( float *x_h, int *nP, float *hist, int *nhistIn, float *xMin, float *xRng, float *weight, float *dx, float *EField )
{
	int nx = *nP;
	float *x_d, *rho_h, *rho_d, *E_h, *E_d;
	unsigned int *hist_h, *hist_d;

	//	check and set device

	#ifndef __DEVICE_EMULATION__
		int GPU_N;
		cudaGetDeviceCount (&GPU_N);
		printf("CUDA-capable device count: %i\n", GPU_N);	
		for (int i=0;i<GPU_N;i++)
		{
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties (&deviceProp,i);
			printf("\nDevice %d: \"%s\"\n", i, deviceProp.name);
		}

		cudaSetDevice ( 0 );
		CUERR
	#endif

	#ifdef SMEM_ATOMICS
		printf("Using shared memory atomics\n");
	#endif
	#ifdef GMEM_ATOMICS
		printf("Using global memory atomics\n");
	#endif

	//	create the possion eqn matrix

	printf ("Creating possion eqn matrix ...\n");
	float *a, *a_d;
	int *ipiv_d;
	a	= (float *)malloc(NHIST*NHIST*sizeof(float));
	//a[0]	= -2;
	//a[1]	=  1;
	//a[NHIST*(NHIST-1)]	=  1;
	//a[NHIST*NHIST-1]	= -2;
	//a[NHIST*NHIST-2]	=  1;
	//for (int i=1;i<NHIST-1;i++){
	//	a[i*NHIST+i-1]	=  1;
	//	a[i*NHIST+1]	= -2;
	//	a[i*NHIST+i+1]	=  1;
	//}
	a[0]	= -2.0;
	a[NHIST]	=  1.0;
	a[NHIST-1]	=  1.0;
	a[NHIST*NHIST-1]	= -2.0;
	a[NHIST*(NHIST-1)-1]	=  1.0;
	for (int i=1;i<NHIST-1;i++){
		a[(i-1)*NHIST+i]	=  1.0;
		a[i*NHIST+i]	= -2.0;
		a[(i+1)*NHIST+i]	=  1.0;
	}
	
	printf ("DONE\n");
	CUERR
	printf ("Allocating a_d on device ...\n");
	cudaMalloc ((void **) &a_d, sizeof(float) * NHIST*NHIST );
	cudaMalloc ((void **) &ipiv_d, sizeof(int) * NHIST );
	CUERR
	printf ("DONE\n");

	cudaMemcpy ( a_d, a, NHIST*NHIST*sizeof(float), cudaMemcpyHostToDevice );
	CUERR

	printf(" nx: %i\n nhist: %i\n xMin: %e\n xRng: %e\n", nx, NHIST, *xMin, *xRng);

	hist_h = (unsigned int *)malloc(sizeof(unsigned int)*NHIST);
	for (int i=0;i<NHIST;i++) hist_h[i] = 0;


	// allocate memory on device

	cudaMalloc ((void **) &x_d, sizeof(float)*nx);
	CUERR
	#ifdef GMEM_ATOMICS
		cudaMalloc ((void **) &hist_d, sizeof(unsigned int)*NHIST);
	#else
		cudaMalloc ((void **) &hist_d, sizeof(unsigned int)*NHISTNBLOCK);
	#endif
	CUERR
	cudaMalloc ((void **) &rho_d, sizeof(float)*NHIST);
	CUERR
	cudaMalloc ((void **) &E_d, sizeof(float)*NHIST);
	CUERR

	// copy from host to device

	cudaMemcpy ( x_d, x_h, sizeof(float)*nx, cudaMemcpyHostToDevice);
	CUERR

	// execute kernel

	printf ("Executing kernels ...\n");
	#ifdef GMEM_ATOMICS
		cudaMemset ( hist_d, 0, sizeof(unsigned int)*NHIST ); // sets all bytes to 0
		hist_kernel<<<NBLOCK,NTHREAD_PB>>>(x_d,nx,hist_d,*xMin,*xRng,rho_d,*weight,*dx);
		remove_average<<<1,1>>>(rho_d,*dx);
	#else
		hist_kernel<<<NBLOCK,NTHREAD_PB>>>(x_d,nx,hist_d,*xMin,*xRng,rho_d,*weight,*dx);
		mergePerBlockHistograms_kernel<<<NHIST,NBLOCK>>>(hist_d,*dx);
		remove_average<<<1,1>>>(rho_d);
	#endif
	printf ("DONE\n");

	rho_h = (float *)malloc(sizeof(float)*NHIST);
	cudaMemcpy ( rho_h, rho_d, sizeof(float)*NHIST, cudaMemcpyDeviceToHost);
	CUERR
	for (int i=0;i<NHIST;i++){
 		 printf ("%i: %f, %f\n",i, rho_h[i] );
	}


	//	initialise cula

	printf ("Initializing cula ...\n");
	culaStatus status;
	status	= culaInitialize();
	checkStatus(status);
	printf ("DONE\n");

	//status	= culaDeviceSgesv ( NHIST, 1, a_d, NHIST, ipiv_d, rho_d, NHIST ); 
	printf ("cula solve status: %i\n", status);
	checkStatus(status);

	//compute_eField<<<1,1>>>(rho_d,E_d,*dx);

	// copy result back to host

	cudaMemcpy ( hist_h, hist_d, sizeof(unsigned int)*NHIST, cudaMemcpyDeviceToHost);
	CUERR

	rho_h = (float *)malloc(sizeof(float)*NHIST);
	cudaMemcpy ( rho_h, rho_d, sizeof(float)*NHIST, cudaMemcpyDeviceToHost);
	CUERR

	E_h = (float *)malloc(sizeof(float)*NHIST);
	cudaMemcpy ( E_h, E_d, sizeof(float)*NHIST, cudaMemcpyDeviceToHost);
	CUERR

	for (int i=0;i<NHIST;i++){
		hist[i] = rho_h[i];
		EField[i]	= E_h[i];	
 		printf ("%i: %f, %f\n",i, rho_h[i], EField[i] );
	}

	printf ("Shutting down cula ...\n");
	culaShutdown();	
	printf ("DONE\n");

	cudaFree ( hist_d );
	cudaFree ( x_d );
	cudaFree ( a_d );
	cudaFree ( E_d );
	CUERR

	//free ( hist_h );
	//free ( E_h );

	return;	
}
