#include "scan2.h"
#include "gpuerrors.h"
#include "stdio.h"
#include "math.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z
__global__ void Kernel1(int *input, int x)
{
    int tid = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;
//int tid = (by * gridDim.x + bx) * blockDim.x  + tx;
	input[tid]=input[tid]+x;



}
__global__ void Kernel(int *input, int *output,int* results, int space, int step, int steps, bool Direction)
{
//int i = (by * gridDim.x + bx) * blockDim.x  + tx;
    int tid = (bz*gridDim.y*gridDim.x + by * gridDim.x + bx) * blockDim.x  + tx;

//int tix = threadIdx.x + blockDim.x * blockIdx.x; 
//int tiy = threadIdx.y + blockDim.y * blockIdx.y;
//int tid = tix + tiy*gridDim.x*blockDim.x;




int res ;
	if(Direction)
	{
		if(tid<space)
		{
			res  = output[tid]; //ONLY REWRITE TO CORRECT MEMORY ADDRESS
			input[tid] = res;
		}
		else
		{
			res = output[tid] + output[tid-space];
			input[tid] = res;
		}
	}
	else
	{
		if(tid<space)
		{
			res = input[tid]; //ONLY REWRITE TO CORRECT MEMORY ADRESS
			output[tid] = res;
		}
		else
		{
			res = input[tid] + input[tid-space];
			output[tid] = res;
		}
	}

if (step == 0) {
results[tid] = -1 * input[tid];
}
if(step == (steps-1)){
results[tid] += res;
}

}


void gpuKernel(int* a, int* c, int n) {	
        int* input1;
	int* output1;
	int* result1;



       
if(n<67108864)
{
	HANDLE_ERROR(cudaMalloc((void**)&input1, n * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&output1, n * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&result1, n * sizeof(int)));
        HANDLE_ERROR(cudaMemcpy(input1, a, n * sizeof(int), cudaMemcpyHostToDevice));

dim3 THREADS_PER_BLOCK(1024,1,1);
	dim3 BLOCKS_PER_GRID( n / 1048576, 32,32);	


//LAUNCH KERNELS IN LOOP
int space = 1;

bool direction=0;
	int steps = (int)log2((float)n);
for ( int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}


	

	HANDLE_ERROR(cudaMemcpy(c, result1, n * sizeof(int), cudaMemcpyDeviceToHost));

        HANDLE_ERROR(cudaFree(input1));
        HANDLE_ERROR(cudaFree(output1));
        HANDLE_ERROR(cudaFree(result1));
}


if(n==67108864)
{
	HANDLE_ERROR(cudaMalloc((void**)&input1, n * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&output1, n * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&result1, n * sizeof(int)));
        HANDLE_ERROR(cudaMemcpy(input1, a, n * sizeof(int), cudaMemcpyHostToDevice));

dim3 THREADS_PER_BLOCK(1024/8,1,1);
	dim3 BLOCKS_PER_GRID( n / 1048576, 32*2,32*4);	




//LAUNCH KERNELS IN LOOP
int space = 1;

bool direction=0;
	int steps = (int)log2((float)n);
for ( int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}


	

	HANDLE_ERROR(cudaMemcpy(c, result1, n * sizeof(int), cudaMemcpyDeviceToHost));

        HANDLE_ERROR(cudaFree(input1));
        HANDLE_ERROR(cudaFree(output1));
        HANDLE_ERROR(cudaFree(result1));
}






if(n==67108864*2)
{
	HANDLE_ERROR(cudaMalloc((void**)&input1,  (n/2) * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&output1, (n/2) * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&result1, (n/2) * sizeof(int)));
	//int* t1;
	//int* c1;
        //t1 = (int*)malloc((n/2) * sizeof(int));
       // c1 = (int*)malloc((n/2) * sizeof(int));
	

//for(int j =0; j< n/2 ;j++)
//	t1[j]=a[j];



        HANDLE_ERROR(cudaMemcpy(input1, a, (n/2) * sizeof(int), cudaMemcpyHostToDevice));


	//dim3 THREADS_PER_BLOCK(1024/8,1,1);
	//dim3 BLOCKS_PER_GRID( n / 1048576, 16,16);	
dim3 THREADS_PER_BLOCK(1024/8,1,1);
	dim3 BLOCKS_PER_GRID( n / (1048576*2), 64,64*2);	

int space = 1;

bool direction=0;
	int steps = (int)log2((float)(n/2));
for ( int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}




//printf("11111");
	HANDLE_ERROR(cudaMemcpy(c, result1, (n/2) * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/2 ;j++)
//	c[j]=t1[j];

int carry = c[(n/2)-1]+a[(n/2)-1];
//////////////////////////////////////////////////
//for(int j =0; j< n/2 ;j++)
//	t1[j]=a[j+(n/2)];

HANDLE_ERROR(cudaMemcpy(input1, &a[n/2], n/2 * sizeof(int), cudaMemcpyHostToDevice));


	//dim3 THREADS_PER_BLOC1(1024,1,1);
	//dim3 BLOCKS_PER_GRID1( n / 1048576,32,32);
 space = 1;

 direction=0;
	 steps = (int)log2((float)n/2);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
//t1=&result1[n/2];
dim3 b(1024,1,1);
	dim3 a( n / (1024*1024),256,1);
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
//Kernel1<<<a,b >>>(result1,carry);
//printf("fghdghdhdhdh");
	HANDLE_ERROR(cudaMemcpy(&c[n/2], result1, n/2 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/2 ;j++)
//	c[j+n/2]+=carry;



        HANDLE_ERROR(cudaFree(input1));
        HANDLE_ERROR(cudaFree(output1));
        HANDLE_ERROR(cudaFree(result1));

}

if(n==67108864*4)
{
	HANDLE_ERROR(cudaMalloc((void**)&input1,  (n/4) * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&output1, (n/4) * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&result1, (n/4) * sizeof(int)));
	
        HANDLE_ERROR(cudaMemcpy(input1, a, (n/4) * sizeof(int), cudaMemcpyHostToDevice));

	
dim3 THREADS_PER_BLOCK(1024/8,1,1);
	dim3 BLOCKS_PER_GRID( n / (1048576*4), 64,64*2);	

int space = 1;

bool direction=0;
	int steps = (int)log2((float)(n/4));
for ( int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}

	HANDLE_ERROR(cudaMemcpy(c, result1, (n/4) * sizeof(int), cudaMemcpyDeviceToHost));

int carry = c[(n/4)-1]+a[(n/4)-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[n/4], n/4 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/4);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);

	HANDLE_ERROR(cudaMemcpy(&c[n/4], result1, n/4 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/4 ;j++)
//	c[j+n/4]+=carry;

 carry = c[(n/2)-1]+a[(n/2)-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[n/2], n/4 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/4);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);

	HANDLE_ERROR(cudaMemcpy(&c[n/2], result1, n/4 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/4 ;j++)
//	c[j+n/2]+=carry;

carry = c[(3*n/4)-1]+a[(3*n/4)-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[3*n/4], n/4 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/4);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);

	HANDLE_ERROR(cudaMemcpy(&c[3*n/4], result1, n/4 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/4 ;j++)
	//c[j+3*n/4]+=carry;



        HANDLE_ERROR(cudaFree(input1));
        HANDLE_ERROR(cudaFree(output1));
        HANDLE_ERROR(cudaFree(result1));

}


if(n==67108864*8)
{

	HANDLE_ERROR(cudaMalloc((void**)&input1,  (n/8) * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&output1, (n/8) * sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&result1, (n/8) * sizeof(int)));
	
        HANDLE_ERROR(cudaMemcpy(input1, a, (n/8) * sizeof(int), cudaMemcpyHostToDevice));

	
dim3 THREADS_PER_BLOCK(1024/8,1,1);
	dim3 BLOCKS_PER_GRID( n / (1048576*8), 64,64*2);	

int space = 1;

bool direction=0;
	int steps = (int)log2((float)(n/8));
for ( int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}

	HANDLE_ERROR(cudaMemcpy(c, result1, (n/8) * sizeof(int), cudaMemcpyDeviceToHost));

int carry = c[(n/8)-1]+a[(n/8)-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[n/8], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
	HANDLE_ERROR(cudaMemcpy(&c[n/8], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+n/8]+=carry;

 carry = c[(n/4)-1]+a[(n/4)-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[n/4], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
	HANDLE_ERROR(cudaMemcpy(&c[n/4], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+n/4]+=carry;

carry = c[(3*n/8)-1]+a[(3*n/8)-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[3*n/8], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);

	HANDLE_ERROR(cudaMemcpy(&c[3*n/8], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+3*n/8]+=carry;

carry = c[(n/8)*4-1]+a[(n/8)*4-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[(n/8)*4], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
	HANDLE_ERROR(cudaMemcpy(&c[(n/8)*4], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+(n/8)*4]+=carry;

carry = c[(n/8)*5-1]+a[(n/8)*5-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[(n/8)*5], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
	HANDLE_ERROR(cudaMemcpy(&c[(n/8)*5], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+(n/8)*5]+=carry;

carry = c[(n/8)*6-1]+a[(n/8)*6-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[(n/8)*6], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
	HANDLE_ERROR(cudaMemcpy(&c[(n/8)*6], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+(n/8)*6]+=carry;

carry = c[(n/8)*7-1]+a[(n/8)*7-1];
//////////////////////////////////////////////////

HANDLE_ERROR(cudaMemcpy(input1, &a[(n/8)*7], n/8 * sizeof(int), cudaMemcpyHostToDevice));

 space = 1;
 direction=0;
 steps = (int)log2((float)n/8);
for (int step = 0; step<steps; step++){
 direction = ((step % 2) !=0);  
Kernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(input1, output1, result1, space,  step, steps,direction);
space = space * 2;
}
	Kernel1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(result1,carry);
	HANDLE_ERROR(cudaMemcpy(&c[(n/8)*7], result1, n/8 * sizeof(int), cudaMemcpyDeviceToHost));

//for(int j =0; j< n/8 ;j++)
//	c[j+(n/8)*7]+=carry;



        HANDLE_ERROR(cudaFree(input1));
        HANDLE_ERROR(cudaFree(output1));
        HANDLE_ERROR(cudaFree(result1));

}







}

