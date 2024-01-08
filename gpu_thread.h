#include "cuda_runtime.h"
#define NUM_OF_THREADS 1024
__global__ void DC(int input_row,
    int input_col,
    int* input,
    int kernel_row,
    int kernel_col,
    int* kernel,
    int output_row,
    int output_col,
    long long unsigned int* output)


{
    int p = (blockIdx.x * blockDim.x) + threadIdx.x;
    int N = output_row * output_col;

    if (p < N) {

        int output_i = p / output_col;
        int output_j = p % output_col;
        for (int kernel_i = 0; kernel_i < kernel_row; kernel_i++)
        {
            for (int kernel_j = 0; kernel_j < kernel_col; kernel_j++)
            {
                int input_i = (output_i + 2 * kernel_i) % input_row;
                int input_j = (output_j + 2 * kernel_j) % input_col;
                output[output_i * output_col + output_j] += input[input_i * input_col + input_j]
                    * kernel[kernel_i * kernel_col + kernel_j];
            }
        }

    }
}
void gpuThread(int input_row,
    int input_col,
    int* input,
    int kernel_row,
    int kernel_col,
    int* kernel,
    int output_row,
    int output_col,
    long long unsigned int* output)
{

    int* cu_in, * cu_ker;
    long long unsigned int* cu_out;
    int N = output_col * output_row;

    cudaMalloc(&cu_in, input_col * input_row * sizeof(int));
    cudaMalloc(&cu_out, output_col * output_row * sizeof(long long unsigned int));
    cudaMalloc(&cu_ker, kernel_col * kernel_row * sizeof(int));

    cudaMemcpy(cu_in, input, input_col * input_row * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(cu_ker, kernel, kernel_col * kernel_row * sizeof(int), cudaMemcpyHostToDevice);

    int N_T = NUM_OF_THREADS;
    int N_B = (N + N_T - 1) / N_T;  

    DC <<<dim3(N_B,1), dim3(N_T,1) >>> (input_row, input_col, cu_in, kernel_row, kernel_col, cu_ker, output_row, output_col, cu_out);

    cudaMemcpy(output, cu_out,  output_col * output_row * sizeof(long long unsigned int), cudaMemcpyDeviceToHost);

    cudaFree(cu_in);
    cudaFree(cu_ker);
    cudaFree(cu_out);

}

