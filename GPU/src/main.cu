#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

// setup settings
__device__ __constant__ size_t d_n = 10; // Hits per Site

// min 2 and has to be the same
__device__ __constant__ size_t d_xjPerThread = 1 << 1;
size_t h_xjPerThread = 1 << 1;

__device__ __constant__ double d_a = 1.0;
__device__ __constant__ double d_lambda = 0.0;
__device__ __constant__ double d_mu_sq = 1.0;
__device__ __constant__ double d_f_sq = 2.0;
__device__ __constant__ double d_m0 = 1.0;
__device__ long seed = 0;
size_t runs = 4; // Runs (increase to yield more data)
size_t SII = 5;  // Statistical Independent Iterations
size_t RS = 0;   // Recording Start
char DATA_PATH[] = "data/data_fig4.csv";

// Guide for Device Setting:
// RTX 2060 super has 34 SMs
// 1 Warp are 32 Threads, 64 Warps per SM at maximum
// So 32*64*34=69,632 Threads can be active at once
// 512 or 1024 Threads per Block are ideal, though at 1024 it is harder to start

// Device dependent settings: do not change
//size_t MCS = 34 * 4 * 1 << 0; 		// 34SMs*4BlocksPerSM=136Blocks: Monte Carlo Simulations => 69,632 Threads => 1 Full GPU
size_t MCS = 1;
size_t MCI = 500;			// 500MCI/5SII=100 Configurations
size_t threadsPerBlock = 1 << 2; 	// 512 threads per block => 4 Blocks = 2048 Threads = 1 Full SM

size_t blocksPerGrid = MCS;

// to increase N, use h_xjPerThread
size_t N = h_xjPerThread * (threadsPerBlock) + 2;

__device__ size_t d_N = (1 << 3) * (1 << 9) + 2;

size_t length = MCS * N;
size_t size = length * sizeof(double);
size_t total_configurations = (MCI - RS) * MCS / SII;

#define gpuErrchk(ans)                                                         \
{ gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
		bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
				line);
		if (abort)
			exit(code);
	}
}
__global__ void init_rand(curandState *state) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	curand_init(42, idx, 0, &state[idx]);
}

__device__ double rand_x(curandState *state, double min, double max) {
	return curand_uniform(state) * (max - min) + min;
}

__device__ double potential_1(double *x) {
	return (d_mu_sq * *x * *x * 0.5 + d_lambda * *x * *x * *x * *x);
}

__device__ double potential_2(double *x) {
	return d_lambda * pow(*x * *x - d_f_sq, 2);
}

__device__ double calc_S_of_xj(double x, double *ptr) {
	return d_m0 * 0.5 * (pow(x - *(ptr + 1), 2) + pow(*(ptr - 1) - x, 2)) / d_a +
		d_a * (potential_1(&x) + potential_1(ptr - 1));
}

__device__ double calc_tot_S(double *sites) {
	double temp = 0;
	for (size_t i = 0; i < d_N / 2 - 1; i++) {
		temp += calc_S_of_xj(sites[i * 2 + 1], sites + i * 2);
	}
	return temp;
}

__device__ double calc_dS(double *xptr, double *new_xptr) {
	return calc_S_of_xj(*new_xptr, xptr) - calc_S_of_xj(*xptr, xptr);
}

__device__ void step(double *xptr, curandState *local_state) {
	double delta = 2 * sqrt(d_a);
	double new_xptr = rand_x(local_state, *xptr - delta, *xptr + delta);
	double dS = calc_dS(xptr, &new_xptr);
	// double r = rand_x(local_state, 0, 1);
	// double edS = pow(M_E, -dS);
	if (dS < 0 || pow(M_E, -dS) > rand_x(local_state, 0, 1))
		*xptr = new_xptr;
}

// cuda c kernel
__global__ void Simulate(double *sites, int iterations, curandState *state) {
	printf("Kernal launched");
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	for (int _ = 0; _ < iterations; _++) {
		for (int __ = 0; __ < d_n; __++) {
			for (int i = 0; i < d_xjPerThread; i++) {
				int site = idx * d_xjPerThread + i + 1;
				step(sites + site, state + idx);
			}
		}
		__syncthreads();
		//if (idx == 0)
		//printf("S: %lf\n", calc_tot_S(sites));
	}
}

__global__ void print_S(double *sites) {}

__global__ void print_sites(double *sites) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	printf("D: site %i, value %lf\n", idx, sites[idx]);
}

void messure(FILE *file, double *values) {
	for (size_t i = 0; i < MCS; i++) {
		for (size_t k = 0; k < N - 1; k++){
			fprintf(file, "%lf;", values[i * N + k]);
		}
		fprintf(file, "%lf\n", values[i * N + N - 1]);
	}
}

int main() {
	printf("Degrees of freedom N: %ld \n", N);

	printf("Launching.\n"); // setup main configs
	FILE *data_file;
	data_file = fopen(DATA_PATH, "w");
	srand(42);
	// size_t set_limit = size_t(2084) * size_t(2084) * size_t(2084) * 2;
	// printf("test %ld", set_limit);
	// size_t limit = 0;
	// cudaDeviceSetLimit(cudaLimitMallocHeapSize, set_limit);
	// cudaDeviceGetLimit(&limit, cudaLimitMallocHeapSize);
	// printf("Availible memory: %ld\n", limit);
	//  const size_t malloc_limit = size_t()

	// setup cuRAND state for randomization
	printf("Setting cuRAND..\n");
	curandState *d_state;
	gpuErrchk(
			cudaMalloc((void **)&d_state, threadsPerBlock * blocksPerGrid * sizeof(curandState)));
	init_rand<<<blocksPerGrid, threadsPerBlock>>>(d_state);
	gpuErrchk(cudaDeviceSynchronize());

	// setup kernel configs
	printf("Threads per block: %ld; Blocks per grid: %ld\n", threadsPerBlock,
			blocksPerGrid);

	// initialise device sites
	double *h_sites, *d_sites;
	gpuErrchk(cudaMallocHost((void **)&h_sites, size));
	gpuErrchk(cudaMalloc((void **)&d_sites, size));
	//h_sites = (double *)malloc(size);


	clock_t start = clock();
	for (int k = 0; k < runs; k++) {
		// initialise host sites and rand
		for (int i = 0; i < length; i++) {
			*(h_sites + i) = ((double)rand() / RAND_MAX - 0.5) * 20;
		}

		// copy sites from host to device
		gpuErrchk(cudaMemcpy(d_sites, h_sites, size, cudaMemcpyHostToDevice));

		// wait until the Simulation comes to equilibrium and take first messurement
		Simulate<<<blocksPerGrid, threadsPerBlock>>>(d_sites, RS, d_state);
		//printf("cuda c err: %s\n", cudaGetLastError());
		cudaGetLastError();
		gpuErrchk(cudaDeviceSynchronize());
		gpuErrchk(cudaMemcpy(h_sites, d_sites, size, cudaMemcpyDeviceToHost));
		messure(data_file, h_sites);

		size_t iterations = MCI / SII - 1;
		for (size_t i = 0; i < iterations; i++) {
			Simulate<<<blocksPerGrid, threadsPerBlock>>>(d_sites, SII, d_state);
			gpuErrchk(cudaDeviceSynchronize());
			gpuErrchk(cudaMemcpy(h_sites, d_sites, size, cudaMemcpyDeviceToHost));
			messure(data_file, h_sites);
			printf("(%3.2lf %%)\n",
					(i + 1 + k * iterations) * 100.0 / iterations / runs);
		}
	}
	clock_t end = clock();
	fclose(data_file);

	// free memory
	printf("Free memory..\n");
	//free(h_sites);
	gpuErrchk(cudaFreeHost(h_sites));
	gpuErrchk(cudaFree(d_sites));

	printf("Done.\n");
	double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
	double benchmark = MCI * MCS * N * runs / time_spend;
	printf("Configurations Messured: %ld; Time Spend on Routine: %0.2lf; "
			"Benchmark: %0.2e sites/s\n",
			total_configurations, time_spend, benchmark);
	return 0;
	;
}
