#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

// setup settings
__device__ __constant__ int d_n = 10; // Hits per Site
__device__ __constant__ double d_a = 1.0;
__device__ __constant__ double d_lambda = 0.0;
__device__ __constant__ double d_mu_sq = 1.0;
__device__ __constant__ double d_f_sq = 2.0;
__device__ __constant__ double d_m0 = 1.0;
int runs = 1 << 2;	// Runs (increase to yield more data)
int SII = 5;    	// Statistical Independent Iterations
int RS = 10;    	// Recording Start

// Device dependent settings: do not change
int MCS = 1 << 10;	// 1024 Monte Carlo Simulations
int MCI = 500;		// 500 Monte Carlo Iterations
int N = 1 << 9;		// 512 degrees of freedom

int size = MCS * N;

__global__ void init_rand(curandState *state) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  curand_init(1234, idx, 0, &state[idx]);
  printf("cuRAND initialized.\n");
}

__device__ double rand_x(curandState *seed, double min, double max) {
  return curand_uniform(seed) * (max - min) + min;
}

__device__ double potential_1(double *x) {
  return (d_mu_sq * *x * *x * 0.5 + d_lambda * *x * *x * *x * *x);
}

__device__ double potential_2(double *x) {
  return d_lambda * pow(*x * *x - d_f_sq, 2);
}

__device__ double calc_S_of_xj(double x, double *ptr) {
  // return d_m0 * 0.5 * pow(x - *(ptr + 1), 2) / d_a + d_a * potential_1(&x) +
  // d_m0 * 0.5 * pow(*(ptr - 1) - x, 2) / d_a + d_a * potential_1(ptr - 1);
  return d_m0 * 0.5 * (pow(x - *(ptr + 1), 2) + pow(*(ptr - 1) - x, 2)) / d_a +
         d_a * (potential_1(&x) + potential_1(ptr - 1));
}

__device__ double calc_dS(double *xptr, double *new_xptr) {
  return calc_S_of_xj(*new_xptr, xptr) - calc_S_of_xj(*xptr, xptr);
}

__device__ void step(double *xptr, curandState *seed) {
  double delta = 2 * sqrt(d_a);
  double new_xptr = rand_x(seed, *xptr - delta, *xptr + delta);
  double dS = calc_dS(xptr, &new_xptr);
  if (dS < 0 || pow(M_E, -dS) > rand_x(seed, 0, 1))
    *xptr = new_xptr;
}

// cuda c kernel
__global__ void Simulate(double *sites, double iterations, curandState *state) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  curandState *seed = state + idx;
  if (threadIdx.x != 0 && threadIdx.x != blockDim.x - 1)
    for (int _ = 0; _ < iterations; _++) {
      // if (idx == 1)
      //  printf("Iteration: %d\n", _);
      for (int __ = 0; __ < d_n; __++) {
        step(sites + idx * 2, seed);
        step(sites + idx * 2 + 1, seed);
      }
      __syncthreads();
    }
}

void messure(FILE *file, double *values) {
  for (size_t i = 0; i < MCS; i++) {
    for (size_t k = 0; k < N - 1; k++)
      fprintf(file, "%lf;", values[i * N + k]);
    fprintf(file, "%lf\n", values[i * N + N - 1]);
  }
}

int main() {

  printf("Launching.\n");
  // setup main configs
  char DATA_PATH[] = "data_fig4.csv";
  FILE *data_file;
  data_file = fopen(DATA_PATH, "w");
  srand(42);

  // setup cuRAND state for randomization
  printf("Setting cuRAND..\n");
  curandState *d_state;
  cudaMalloc(&d_state, sizeof(curandState));
  init_rand<<<1, 1>>>(d_state);
  cudaDeviceSynchronize();

  // setup kernel configs
  int threadsPerBlock = int(N / 2);
  int blocksPerGrid = MCS;
  printf("Threads per block: %d; Blocks per grid: %d\n", threadsPerBlock,
         blocksPerGrid);
  cudaDeviceSynchronize();

  // initialise device sites
  double *d_sites;
  cudaMalloc(&d_sites, size);
  double *h_sites;

  clock_t start = clock();
  for (int k = 0; k < runs; k++) {
    // initialise host sites and rand
    printf("Initialising host sites..\n");
    h_sites = (double *)malloc(size * sizeof(double));
    for (size_t i = 0; i < size; i++) {
      h_sites[i] = ((double)rand() / RAND_MAX - 0.5) * 2;
    }

    // copy sites from host to device
    cudaMemcpy(d_sites, h_sites, size, cudaMemcpyHostToDevice);

    // wait until the Simulation comes to equilibrium and take first messurement
    printf("Reaching equilibrium..\n");
    Simulate<<<MCS, threadsPerBlock>>>(d_sites, RS, d_state);
    cudaDeviceSynchronize();
    cudaMemcpy(h_sites, d_sites, size, cudaMemcpyDeviceToHost);
    messure(data_file, h_sites);

    printf("Messuring..\n");
    // Simulate and take messurements
    size_t iterations = MCI / SII - 1;
    for (size_t i = 0; i < iterations; i++) {
      Simulate<<<blocksPerGrid, threadsPerBlock>>>(d_sites, SII, d_state);
      cudaDeviceSynchronize();

      cudaMemcpy(h_sites, d_sites, size, cudaMemcpyDeviceToHost);
      messure(data_file, h_sites);

      printf("(%3.2lf %%)\n", (i + 1) * 100.0 / iterations);
    }
  }
  clock_t end = clock();
  fclose(data_file);

  // free memory
  printf("Free memory..\n");
  free(h_sites);
  cudaFree(d_sites);

  printf("Done.\n");
  double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
  double benchmark = MCI * MCS * N * runs / time_spend;
  printf("Time Spend on Routin: %lf; Benchmark: %e sites/s\n", time_spend,
         benchmark);
  return 0;
  ;
}
