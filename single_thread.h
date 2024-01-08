#include<stdlib.h>
#include<stdio.h>
#include<immintrin.h>
//optimize this function
//1) used avx instructions and other optimizations.

void singleThread(int input_row,
  int input_col,
  int * input,
  int kernel_row,
  int kernel_col,
  int * kernel,
  int output_row,
  int output_col,
  long long unsigned int * output) 



/*{

    for(int i = 0; i < output_row * output_col; ++i)
        output[i] = 0;

    for(int output_i = 0; output_i< output_row; output_i++)
    {
        for(int output_j = 0; output_j< output_col; output_j++)
        {
            for(int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
            {
                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j++)
                {
                    int input_i = (output_i + 2*kernel_i) % input_row;
                    int input_j = (output_j + 2*kernel_j) % input_col;
                    output[output_i * output_col + output_j] += input[input_i*input_col +input_j] * kernel[kernel_i*kernel_col +kernel_j];
                }
            }
        }
    }

}*/

{
//get the no of rows and cols that can be calculated contiguously
    int remain_rows=output_row-kernel_row+1;
    int remain_cols=output_col-kernel_col+1;

    remain_rows=remain_rows-(remain_rows%8);
    remain_cols=remain_cols-(remain_cols%8);

    int * t=(int *)malloc(sizeof(int)*output_row*output_col);
    for(int kernel_i=0;kernel_i<kernel_row;kernel_i++)
    {
        for(int kernel_j=0;kernel_j<kernel_col;kernel_j++)
        {
            __m256i kernel_values=_mm256_set1_epi32(kernel[kernel_i*kernel_col+kernel_j]);
            for(int output_i=0;output_i<output_row;output_i++)
            {
                int input_i=(output_i+2*kernel_i)%input_row;
                int inp_idx=input_i*input_col;
                int otp_idx=output_i*output_col;

                for(int output_j=0;output_j<remain_cols;output_j+=8)
                {
                    int input_j=(output_j+2*kernel_j);
                    __m256i input_values = _mm256_loadu_si256((__m256i *)(input+inp_idx+input_j));
                    __m256i result_values = _mm256_mullo_epi32(kernel_values,input_values);

                    result_values=_mm256_add_epi32(result_values,_mm256_loadu_si256((__m256i *)(t+output_i*output_col+output_j)));
                    _mm256_storeu_si256((__m256i *)(t+output_i*output_col+output_j),result_values);
                }
            }
        }
    }


    for(int i = 0; i < output_row * output_col; ++i)
    output[i] = t[i];

    /*for(int output_i = 0; output_i< output_row; output_i++)
    {
        for(int output_j = remain_cols; output_j< output_col; output_j++)
        {
            for(int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
            {
                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j++)
                {
                    int input_i = (output_i + 2*kernel_i) % input_row;
                    int input_j = (output_j + 2*kernel_j) % input_col;
                    output[output_i * output_col + output_j] += input[input_i*input_col +input_j] * kernel[kernel_i*kernel_col +kernel_j];
                }
            }
        }
    }*/

    for(int output_i = 0; output_i< output_row; output_i++)
    {
        int a = output_i * output_col;
        for(int output_j = remain_cols; output_j< output_col; output_j++)
        {
            int out_idx = a + output_j;
            long long unsigned int x = 0;
            for(int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
            {
                int input_i = output_i + (kernel_i)*2;
                if(input_i >= input_row)
                {
                    input_i -= input_row;
                }
                int b = input_i*input_col;
                int c = kernel_i*kernel_col;

                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j= kernel_j+4)
                {
                    int input_j = output_j + (kernel_j)*2;
                    if(input_j >= input_col)
                    {
                        input_j -= input_col;
                    }
                    x += input[b +input_j] * kernel[c + kernel_j];
                    if(kernel_j + 1 < kernel_col)
                    {
                        input_j = output_j + 2*(kernel_j+1);
                        if(input_j >= input_col)
                        {
                            input_j -= input_col;
                        }
                        x += input[b + input_j] * kernel[c + kernel_j+1];
                    }
                    if(kernel_j + 2 < kernel_col)
                    {
                        input_j = output_j + 2*(kernel_j+2);
                        if(input_j >= input_col)
                        {
                            input_j -= input_col;
                        }
                        x += input[b + input_j] * kernel[c + kernel_j+2];
                    }
                    if(kernel_j + 3 < kernel_col)
                    {
                        input_j = output_j + 2*(kernel_j+3);
                        if(input_j >= input_col)
                        {
                            input_j -= input_col;
                        }
                        x += input[b + input_j] * kernel[c + kernel_j+3];
                    }

                }
            }
            output[out_idx] = x;
        }
    }

   
}
/*{

    for(int output_i = 0; output_i< output_row; output_i++)
    {
        int a = output_i * output_col;
        for(int output_j = 0; output_j< output_col; output_j++)
        {
            int out_idx = a + output_j;
            long long unsigned int x = 0;
            for(int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
            {
                int input_i = output_i + (kernel_i)*2;
                if(input_i >= input_row)
                {
                    input_i -= input_row;
                }
                int b = input_i*input_col;
                int c = kernel_i*kernel_col;

                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j= kernel_j+4)
                {
                    int input_j = output_j + (kernel_j)*2;
                    if(input_j >= input_col)
                    {
                        input_j -= input_col;
                    }
                    x += input[b +input_j] * kernel[c + kernel_j];
                    if(kernel_j + 1 < kernel_col)
                    {
                        input_j = output_j + 2*(kernel_j+1);
                        if(input_j >= input_col)
                        {
                            input_j -= input_col;
                        }
                        x += input[b + input_j] * kernel[c + kernel_j+1];
                    }
                    if(kernel_j + 2 < kernel_col)
                    {
                        input_j = output_j + 2*(kernel_j+2);
                        if(input_j >= input_col)
                        {
                            input_j -= input_col;
                        }
                        x += input[b + input_j] * kernel[c + kernel_j+2];
                    }
                    if(kernel_j + 3 < kernel_col)
                    {
                        input_j = output_j + 2*(kernel_j+3);
                        if(input_j >= input_col)
                        {
                            input_j -= input_col;
                        }
                        x += input[b + input_j] * kernel[c + kernel_j+3];
                    }

                }
            }
            output[out_idx] = x;
        }
    }
}*/
/*{
  {
    for (int i = 0; i < output_row * output_col; ++i)
      output[i] = 0;

    int * temp_out = (int * ) malloc(sizeof(int) * output_row * output_col); // temparory output for output matrix

    int con_rows = output_row - kernel_col + 1 - ((output_row - kernel_row + 1) % 8); // scaling rows of output to multiple of 8
    int con_cols = output_col - kernel_col + 1 - ((output_col - kernel_col + 1) % 8); // scaling cols of output to multiple of 8

    for (int kernel_i = 0; kernel_i < kernel_row; kernel_i++) {
      int kernel__ = kernel_i * kernel_col;
      for (int kernel_j = 0; kernel_j < kernel_col; kernel_j++) {

        __m256i ker_vec = _mm256_set1_epi32(kernel[kernel__ + kernel_j]);

        for (int output_i = 0; output_i < output_row; output_i++) {
          int output__ = output_i * output_col;
          int input_i = (output_i + 2 * kernel_i) % input_row;
          int input__ = input_i * input_col;
          for (int output_j = 0; output_j < con_cols; output_j += 8) {
            int input_j = (output_j + 2 * kernel_j); // no need of % because columns wont wrap around here.
            __m256i inp_vec = _mm256_loadu_si256((__m256i * )(input + input__ + input_j));
            __m256i res_vec = _mm256_mullo_epi32(ker_vec, inp_vec);
            res_vec = _mm256_add_epi32(res_vec, _mm256_loadu_si256((__m256i * )(temp_out + output_i * output_col + output_j)));
            _mm256_storeu_si256((__m256i * )(temp_out + output_i * output_col + output_j), res_vec);
          }
        }
      }
    }

    for (int i = 0; i < output_row * output_col; i++) {
      output[i] = temp_out[i]; // copy temp_out to ouput matrix
    }

    for (int output_i = 0; output_i < output_row; output_i++) {
      int output__ = output_i * output_col;
      for (int output_j = con_cols; output_j < output_col; output_j++) // computing output for remaining cells. (remaining cells can be atmost 7 in each row.)
      {
        unsigned long long temp = 0;
        for (int kernel_i = 0; kernel_i < kernel_row; kernel_i++) {
          int kernel__ = kernel_i * kernel_col;
          for (int kernel_j = 0; kernel_j < kernel_col; kernel_j++) {
            int input_i = (output_i + 2 * kernel_i) % input_row;
            int input_j = (output_j + 2 * kernel_j) % input_col;
            temp += input[input_i * input_col + input_j] * kernel[kernel__ + kernel_j];
          }
        }
        output[output__ + output_j] = temp;
      }
    }
  }
}

*/



////2) used loop unrolling and other optimizations.
// //Optimize this function
//void singleThread(int input_row, 
//                int input_col,
//                int *input, 
//                int kernel_row, 
//                int kernel_col, 
//                int *kernel,
//                int output_row, 
//                int output_col, 
//                long long unsigned int *output ) 
//{
    
//    for(int i = 0; i < output_row * output_col; ++i)
//        output[i] = 0;
//
//    for(int output_i = 0; output_i< output_row; output_i++)
//    {
//        int output__ = output_i * output_col;
//        for(int output_j = 0; output_j< output_col; output_j+=2)     //output is unrolled by factor 2
//        { 
//		  long long  unsigned  int temp_out = 0;
//		  long long  unsigned  int temp_out1 = 0;
//		  int out_temp_out=output__ + output_j;
//			if(output_j < output_col-1)
//            {
//            for( int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
//            {   
//		    	int input__=((output_i + (kernel_i<<1)) % input_row)*input_col;
//                int kernel__=kernel_i*kernel_col;
//                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j+=4)     //kernel is unrolled by factor 4
//                {  
//                    int kernel_ind = kernel__ +kernel_j;
//                    int input_j  = (output_j + ((kernel_j)<<1)) % input_col;
//                    int input_j1 = (output_j + ((kernel_j+1)<<1)) % input_col;
//                    int input_j2 = (output_j + ((kernel_j+2)<<1)) % input_col;
//                    int input_j3 = (output_j + ((kernel_j+3)<<1)) % input_col;
//                    if(kernel_j < kernel_col-3)
//                    temp_out += input[input__+input_j] * kernel[kernel_ind]
//                                                             +input[input__ +input_j1] * kernel[kernel_ind+1]
//                                                             +input[input__ +input_j2] * kernel[kernel_ind+2]
//                                                             +input[input__ +input_j3] * kernel[kernel_ind+3];
//                   	else if(kernel_j < kernel_col -2 )
//                    temp_out += input[input__ +input_j] * kernel[kernel_ind]
//                                                             +input[input__ +input_j1] * kernel[kernel_ind+1]
//                                                             +input[input__ +input_j2] * kernel[kernel_ind+2];
//                  
//                   	else if(kernel_j < kernel_col -1)
//                    temp_out += input[input__+input_j] * kernel[kernel_ind]
//                                                             +input[input__+input_j1] * kernel[kernel_ind+1];
//                    else
//                    temp_out += input[input__+input_j] * kernel[kernel_ind];                                              
//                }
//             }
//				
//		    for(int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
//            {   
//				int input__=((output_i + (kernel_i<<1)) % input_row)*input_col;
//                int kernel__=kernel_i*kernel_col;
//                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j+=4)
//                {  
//                    int kernel_ind = kernel__ +kernel_j;
//                    int input_j = (output_j+1 + (kernel_j<<1)) % input_col;
//                    int input_j1 = (output_j+1 + ((kernel_j+1)<<1)) % input_col;      
//                    int input_j2 = (output_j+1+ ((kernel_j+2)<<1)) % input_col;
//                    int input_j3 = (output_j+1+ ((kernel_j+3)<<1)) % input_col;
//                    if(kernel_j < kernel_col-3)
//                    temp_out1 += input[input__+input_j] * kernel[kernel_ind]
//                                                             +input[input__ +input_j1] * kernel[kernel_ind+1]
//                                                             +input[input__ +input_j2] * kernel[kernel_ind+2]
//                                                             +input[input__ +input_j3] * kernel[kernel_ind+3];
//
//                   else if(kernel_j < kernel_col -2 )
//                    temp_out1 += input[input__ +input_j] * kernel[kernel_ind]
//                                                             +input[input__ +input_j1] * kernel[kernel_ind+1]
//                                                             +input[input__ +input_j2] * kernel[kernel_ind+2];
//                   else if(kernel_j < kernel_col -1)
//                    temp_out1 += input[input__+input_j] * kernel[kernel_ind]
//                                                             +input[input__+input_j1] * kernel[kernel_ind+1];
//                  else
//                    temp_out1 += input[input__+input_j] * kernel[kernel_ind];                                      
//                }
//            }   
//            output[out_temp_out]=temp_out;
//		    output[out_temp_out+1]=temp_out1;
//			}
//			else
//			{			
//            for( int kernel_i = 0; kernel_i< kernel_row; kernel_i++)
//            {   
//				int input__=((output_i + (kernel_i<<1)) % input_row)*input_col;
//                int kernel__=kernel_i*kernel_col;
//                for(int kernel_j = 0; kernel_j< kernel_col; kernel_j+=4)
//                {  
//                    int kernel_ind = kernel__ +kernel_j;
//					int input_j = (output_j + (kernel_j<<1)) % input_col;
//                    int input_j1 = (output_j + ((kernel_j+1)<<1)) % input_col;
//                    int input_j2 = (output_j + ((kernel_j+2)<<1)) % input_col;
//                    int input_j3 = (output_j + ((kernel_j+3)<<1)) % input_col; 
//                    if(kernel_j < kernel_col-3)
//                    temp_out += input[input__+input_j] * kernel[kernel_ind]
//                                                             +input[input__ +input_j1] * kernel[kernel_ind+1]
//                                                             +input[input__ +input_j2] * kernel[kernel_ind+2]
//                                                             +input[input__ +input_j3] * kernel[kernel_ind+3];
//                   else if(kernel_j < kernel_col -2 )
//                    temp_out += input[input__ +input_j] * kernel[kernel_ind]
//                                                             +input[input__ +input_j1] * kernel[kernel_ind+1]
//                                                             +input[input__ +input_j2] * kernel[kernel_ind+2];
//                   else if(kernel_j < kernel_col -1)
//                    temp_out += input[input__+input_j] * kernel[kernel_ind]
//                                                             +input[input__+input_j1] * kernel[kernel_ind+1];
//
//                  else
//                    temp_out += input[input__+input_j] * kernel[kernel_ind];                                      
//
//                }
//            }
//            output[out_temp_out]=temp_out;
//			}
//        }
//     }
//}

