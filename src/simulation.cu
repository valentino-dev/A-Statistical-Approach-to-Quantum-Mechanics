#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

// Guide for Device Setting:
// RTX 2060 super has 34 SMs
// 1 Warp are 32 Threads, 64 Warps per SM at maximum
// So 32*64*34=69,632 Threads can be active at once
// 512 or 1024 Threads per Block are ideal, though at 1024 it is harder to start

#define gpuErrchk(ans)                                                         \
{ gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort)
			exit(code);
	}
}

typedef struct{
	int runs, SII, RS, xjPerThread, n, threadsPerBlock, MCS, MCI, N, total_configurations, blocksPerGrid, array_size;
	double a, m0, lambda, mu_sq, f_sq, delta;
} setting;


__global__ void init_rand(curandState *state) {
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
	curand_init(42, idx, 0, &state[idx]);
}

__device__ double rand_x(curandState *local_state, double min, double max) {
	// return 0.5;
	return curand_uniform(local_state) * (max - min) + min;
}

__device__ double potential_1(double *x, setting *settings) {
	return settings->mu_sq * pow(*x, 2) * 0.5 + settings->lambda * pow(*x, 4);
}

__device__ double potential_2(double *x, setting *settings) {
	return settings->lambda * pow(*x * *x - settings->f_sq, 2);
}

__device__ double calc_S_of_xj(double x, double *ptr, double *previous_site, double *following_site, setting *settings) {
	return settings->m0 * 0.5 * (pow(*following_site - x, 2) + pow(x - *previous_site, 2)) / settings->a +
		settings->a * (potential_1(&x, settings) + potential_1(ptr - 1, settings));
}

__device__ double calc_dS(double *xptr, double *previous_site, double *following_site, double *new_xptr, setting *settings) {
	return calc_S_of_xj(*new_xptr, xptr, previous_site, following_site, settings) - 
		calc_S_of_xj(*xptr, xptr, previous_site, following_site, settings);
}

__device__ void step(double *sites, int site, int previous_site, int following_site, curandState local_state, setting *settings) {
	double *xptr = sites+site;
	double new_x = rand_x(&local_state, *xptr - settings->delta, *xptr + settings->delta);
	
	//printf("a new_xptr: %lf, xptr: %lf, delta: %lf\n", new_x, *xptr, settings->delta);
	//printf("new_x: %lf, old_x: %lf\n", new_x, *xptr);
	double dS = calc_dS(xptr, sites+previous_site, sites+following_site, &new_x, settings);
	//printf("S: %lf\n", dS);
	if (dS < 0 || pow(M_E, -dS) > curand_uniform(&local_state)){
		*xptr = new_x;
		//printf("accept S\n");
	}
}

__device__ void Print_Action(double *sites, setting *settings){
	double S = 0;
	for (size_t i = 0; i < settings->N; i++)
		S += settings->m0 * 0.5 * (pow(sites[(i+1)%settings->N] - sites[i], 2)) / settings->a + settings->a * potential_1(&sites[i], settings);

	printf("S: %lf\n", S);
	
}

// cuda c kernel
__global__ void Simulate(double *sites, int iterations, curandState *state, setting *settings) {
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

	for (int _ = 0; _ < iterations; _++) {
		for (int i = 0; i < settings->xjPerThread; i++) {
			int site = idx * settings->xjPerThread + i;
			int previous_site = (threadIdx.x * settings->xjPerThread + i - 1) % settings->N + settings->N * blockIdx.x;
			int following_site = (threadIdx.x * settings->xjPerThread + i + 1) % settings->N + settings->N * blockIdx.x;
			for (int __ = 0; __ < settings->n; __++) {
				step(sites, site, previous_site, following_site, state[idx], settings);
			}
		}
		__syncthreads();
	}
}

__global__ void print_sites(double *sites) {
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
	printf("D: site %i, value %lf\n", idx, sites[idx]);
}

void messure(FILE *file, double *values, setting *settings) {
	for (size_t i = 0; i < settings->MCS; i++) {
		for (size_t k = 0; k < settings->N - 1; k++){
			fprintf(file, "%lf;", values[i * settings->N + k]);
		}
		fprintf(file, "%lf\n", values[(i + 1) * settings->N - 1]);
	}
}

__global__ void initial_ensamble(double *sites, curandState *state, setting *settings){
	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
	for (size_t i = 0; i < settings->n; i++) {
		sites[idx * settings->n + i] = rand_x(state + idx, -100, 100);
	}
}

void collectSettings(setting *settings, char **argv){
	settings->runs 				= atof(argv[2]);
	settings->SII 				= atof(argv[3]);
	settings->RS 				= atof(argv[4]);
	settings->xjPerThread 		= atof(argv[5]);
	settings->n 				= atof(argv[6]);
	settings->threadsPerBlock 	= atof(argv[7]);
	settings->MCS 				= atof(argv[8]);
	settings->MCI 				= atof(argv[9]);
	settings->a 				= atof(argv[10]);
	settings->m0 				= atof(argv[11]);
	settings->lambda 			= atof(argv[12]);
	settings->mu_sq 			= atof(argv[13]);
	settings->f_sq 				= atof(argv[14]);

	settings->total_configurations = (int) (settings->MCS * settings->runs * settings->MCI / settings->SII);
	settings->N = settings->xjPerThread * settings->threadsPerBlock;
	settings->array_size = settings->MCS * settings->N * sizeof(double);
	settings->blocksPerGrid = settings->MCS;
	settings->delta = 2*sqrt(settings->a);

}


int main(int argc, char** argv) {
	// Checking arguments for error
	if (argc != 15){
		printf("ERROR: Not enough arguments!\nPlease specify the following settings: data_path, runs, SII, RS, xjPerThread, n, threadsPerBlock, MCS, MCI, a, m0, lambda, mu_sq, f_sq.\n");
		return 0;
	}
	printf("%d Arguments:", argc);
	for(size_t i = 0; i < argc; i++){
		printf(" %i: %s,", i, argv[i]);
	}
	printf("\n\n");


	// Collecting Settings: data_path, runs, SII, RS, xjPerThread, n, threadsPerBlock, MCS, MCI, a, m0, lambda, mu_sq, f_sq
	printf("Collecting Settings..\n"); 
	setting *h_settings, *d_settings;
	gpuErrchk(cudaMallocHost((void **)&h_settings, sizeof(setting)));
	gpuErrchk(cudaMalloc((void **)&d_settings, sizeof(setting)));
	collectSettings(h_settings, argv);
	gpuErrchk(cudaMemcpy(d_settings, h_settings, sizeof(setting), cudaMemcpyHostToDevice));
	printf("Settings: %i runs, %i SII, %i RS, %i xjPerthread, n=%i, %i threadsPerBlock, %i MCS, %i MCI, a=%lf, m0=%lf, lambda=%lf, mu_sq=%lf, f_sq=%lf\n", h_settings->runs, h_settings->SII, h_settings->RS, h_settings->xjPerThread, h_settings->n, h_settings->threadsPerBlock, h_settings->MCS, h_settings->MCI, h_settings->a, h_settings->m0, h_settings->lambda, h_settings->mu_sq, h_settings->f_sq);
	printf("Degrees of freedom N: %ld \n", h_settings->N);
	printf("Threads per block: %d; Blocks per grid: %d\n", h_settings->threadsPerBlock, h_settings->blocksPerGrid);


	// Setup main settings
	printf("Launching.\n"); 
	FILE *data_file;
	data_file = fopen(argv[1], "w");
	srand(42);


	// Setup cuRAND state for randomization
	printf("Setting cuRAND..\n");
	curandState *d_state;
	gpuErrchk(cudaMalloc((void **)&d_state, h_settings->threadsPerBlock * h_settings->blocksPerGrid * sizeof(curandState)));
	init_rand<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_state);
	gpuErrchk(cudaDeviceSynchronize());


	// Initialise device sites
	double *h_sites, *d_sites;
	gpuErrchk(cudaMallocHost((void **)&h_sites, h_settings->array_size));
	gpuErrchk(cudaMalloc((void **)&d_sites, h_settings->array_size));
	gpuErrchk(cudaMemset(h_sites, 0, h_settings->array_size));
	gpuErrchk(cudaMemset(d_sites, 0, h_settings->array_size));


	// Starting runs
	clock_t start = clock();
	for (size_t k = 0; k < h_settings->runs; k++) {
		gpuErrchk(cudaMemcpy(d_sites, h_sites, h_settings->array_size, cudaMemcpyHostToDevice));
		initial_ensamble<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_sites, d_state, d_settings);


		// Wait until the Simulation comes to equilibrium
		Simulate<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_sites, h_settings->RS, d_state, d_settings);


		// Simulating and messuring
		size_t iterations = h_settings->MCI / h_settings->SII;
		for (size_t i = 0; i < iterations; i++) {
			gpuErrchk(cudaDeviceSynchronize());
			gpuErrchk(cudaMemcpy(h_sites, d_sites, h_settings->array_size, cudaMemcpyDeviceToHost));
			Simulate<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_sites, iterations, d_state, d_settings);
			messure(data_file, h_sites, h_settings);
			printf("(%3.2lf %%)\n", (i + 1 + k * iterations) * 100.0 / iterations / h_settings->runs);
		}
	}
	clock_t end = clock();
	fclose(data_file);
	gpuErrchk(cudaDeviceSynchronize());


	// Free memory
	printf("Free memory..\n");
	gpuErrchk(cudaFreeHost(h_sites));
	gpuErrchk(cudaFree(d_sites));


	// End Results
	printf("Done.\n");
	double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
	double benchmark = h_settings->MCI * h_settings->MCS * h_settings->N * h_settings->runs / time_spend;
	printf("Configurations Messured: %ld; Time Spend on Routine: %0.2lfs; "
			"Benchmark: %0.2e sites/s\n",
			h_settings->total_configurations, time_spend, benchmark);
	return 0;
}
