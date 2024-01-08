# two-dilated-convolution
Dilated Convolutions in CNN(Convolutional Neural networks) are a type of convolution that “inflate” the kernel by inserting holes between the kernel elements. An additional parameter (dilation rate) indicates how much the kernel is widened.

In this project I implemented 2 dilated convolution alogorithm and optimized it using various techniques.

1)frequency reduction + loop unrolling + various other.

2)Above all mentioned + multi threading in C++(8 threads used in code).

3)GPU based code 

Results:

Speedup achieved for each implentation is :

1) 4.5 for optimizations
2) 9.0 for multi threading
3) 150.0 for GPU implementation.

for full repository and code visit : https://github.com/Abhishekghosh1998/hpca-assignment-2023

How to run:

1)Clone the above github repository 

2)replace  /PartA/header/single_thread.h , /PartA/header/multi_thread.h  and PartB/header/gpu_thread.h files with above files.

3)running the program and everything can be seen in the readme file of that repository only.

