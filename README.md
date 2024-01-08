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

for full repository and code visit : link
