
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <stdio.h>
#include <time.h>


__global__ void printFromGPU(){
    printf("Hello World from GPU\n");
}


int main(int argc, char **argv){
    printFromGPU<<<1, 1>>>();
    return 0;
}
