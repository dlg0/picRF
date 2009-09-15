#include <stdio.h>
#include "cuda.h"


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

static __global__ void hist_kernel ( float *xIn, int nx, unsigned int *histOutBlock, float xMin, float xRng)
{
	const int myThrdId = blockIdx.x * blockDim.x + threadIdx.x;
	const int nThreads = blockDim.x * gridDim.x;
	unsigned int ii;
	float inv_xRng;

	inv_xRng	= 1.0 / xRng;

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
		ii = lrintf((xIn[i]-xMin) * inv_xRng * (NHIST-1));

		// increment my warp histogram, remember there are 6 per block
		warp_hist_kernel ( s_warp_hist+warpOffset, ii, threadTag );

	}

	#ifdef __DEVICE_EMULATION__
		int sum = 0;
		for (int i=0;i<NHIST;i++){
			sum += s_warp_hist[i+warpOffset];
			printf ("%i: %i, %i\n",i, s_warp_hist[i+warpOffset], sum );
		}
	#endif

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

	#ifdef __DEVICE_EMULATION__
		sum = 0;
		for (int i=0;i<NHIST;i++){
			sum += histOutBlock[i+blockIdx.x*NHIST];
			printf ("%i: %i, %i\n",i, histOutBlock[i+blockIdx.x*NHIST], sum );
		}
	#endif

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

/*
	fortran calls this subroutine with
	arguments passed by reference
*/

extern "C" void cudahist_ ( float *x, int *nP, int *hist, int *nhistIn, float *xMin, float *xRng )
{
	int nx = *nP;
	float *x_h, *x_d;
	unsigned int *hist_h, *hist_d;

	//	check and set device

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

	#ifdef SMEM_ATOMICS
		printf("Using shared memory atomics\n");
	#endif
	#ifdef GMEM_ATOMICS
		printf("Using global memory atomics\n");
	#endif


	x_h = (float *)malloc(sizeof(float)*nx);
	for (int i=0;i<nx;i++) x_h[i] = x[i];

	//for (int i=0;i<nx;i++)
	//	printf("%e %e\n",x[i], x_h[i]);

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
	
	// copy from host to device

	cudaMemcpy ( x_d, x, sizeof(float)*nx, cudaMemcpyHostToDevice);
	CUERR

	// execute kernel

	#ifdef GMEM_ATOMICS
		cudaMemset ( hist_d, 0, sizeof(unsigned int)*NHIST ); // sets all bytes to 0
		hist_kernel<<<NBLOCK,NTHREAD_PB>>>(x_d,nx,hist_d,*xMin,*xRng);
	#else
		hist_kernel<<<NBLOCK,NTHREAD_PB>>>(x_d,nx,hist_d,*xMin,*xRng);
		mergePerBlockHistograms_kernel<<<NHIST,NBLOCK>>>(hist_d);
	#endif

	// copy result back to host

	cudaMemcpy ( hist_h, hist_d, sizeof(unsigned int)*NHIST, cudaMemcpyDeviceToHost);
	CUERR

	cudaFree ( hist_d );
	CUERR
	cudaFree ( x_d );
	CUERR

	for (int i=0;i<NHIST;i++) hist[i] = hist_h[i];

	free ( hist_h );

	//int sum = 0;
	//for (int i=0;i<NHIST;i++){
	//	sum += hist[i];
	//	printf ("%i, %i\n",hist[i], sum );
	//}
	
	return;	
}
