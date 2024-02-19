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
inline void gpuAssert(cudaError_t code, const char *file, int line,
		bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
				line);
		if (abort)
			exit(code);
	}
}

typedef struct{
	int runs, SII, RS; xjPerTHread, n, threadsPerBlock, MCS, MCI, a, m0, lmabda, mu_sq, f_sq, delta, N, total_configurations, blocksPerGrid, array_size;
} setting;


__global__ void init_rand(curandState *state) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	curand_init(42, idx, 0, &state[idx]);
		printf("test\n");
}

__device__ double rand_x(curandState *state, double min, double max) {
	return curand_uniform(state) * (max - min) + min;
}

__device__ double potential_1(double *x, settings *settings) {
	return settings->mu_sq * pow(*x, 2) * 0.5 + settings->lambda * pow(*x, 4);
}

__device__ double potential_2(double *x, settings *settings) {
	return settings->lambda * pow(*x * *x - settings->f_sq, 2);
}

__device__ double calc_S_of_xj(double x, double *ptr, double *previous_site, double *following_site, setting *settings) {
	return settings->m0 * 0.5 * (pow(*following_site - x, 2) + pow(x - *previous_site, 2)) / settings->a +
		settings->a * (potential_1(&x, settings) + potential_1(ptr - 1, settings));
}

__device__ double calc_dS(double *xptr, double *previous_site, double *following_site, double *newx_ptr, setting *settings) {
	return calc_S_of_xj(*newx_ptr, xptr, previous_site, following_site, settings) - calc_S_of_xj(*xptr, xptr, previous_site, following_site, settings);
}

__device__ void step(double *sites, int site, int previous_site, int following_site, curandState *local_state, setting *settings) {
	double *xptr = sites+site;
	double new_x = rand_x(local_state, *xptr - settings->delta, *xptr + settings->delta);
	double dS = calc_dS(xptr, sites+previous_site, sites+following_site, &new_x, settings);
	if (dS < 0 || pow(M_E, -dS) > rand_x(local_state, 0, 1))
		*xptr = new_x;
}

__device__ void Print_Action(double *sites, setting *settings){
	double S = 0;
	for (int i = 0; i < settings->N; i++){	
		S += settings->9 * 0.5 * (pow(sites[(i+1)%N] - sites[i], 2)) / settings->a + settings->a * potential_1(&sites[i], settings);
	}

	printf("S: %lf\n", S);
	
}

// cuda c kernel
__global__ void Simulate(double *sites, int iterations, curandState *state, setting *settings) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	for (int _ = 0; _ < iterations; _++) {
		for (int i = 0; i < settings->d_xjPerThread; i++) {
			int site = idx * settings->d_xjPerThread + i;
			int previous_site = (threadIdx.x * settings->d_xjPerThread + i - 1) % (blockDim.x * settings->d_xjPerThread) + blockDim.x * settings->d_xjPerThread * blockIdx.x;
			int following_site = (threadIdx.x * settings->d_xjPerThread + i + 1) % (blockDim.x * settings->d_xjPerThread) + blockDim.x * settings->d_xjPerThread * blockIdx.x;
			for (int __ = 0; __ < settings[4]; __++) {
				step(sites, site, previous_site, following_site, state + idx, settings);
			}
		}
		__syncthreads();
		printf("test\n");
		//if(idx == 0){
			//Print_Action(sites, settings);
				
		//}
	}
}

__global__ void print_sites(double *sites) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	printf("D: site %i, value %lf\n", idx, sites[idx]);
}

void messure(FILE *file, double *values, setting *settings) {
	for (size_t i = 0; i < settings->MCS; i++) {
		for (size_t k = 0; k < setting->N - 1; k++){
			fprintf(file, "%lf;", values[i * N + k]);
		}
		fprintf(file, "%lf\n", values[(i + 1) * N - 1]);
	}
}

__global__ void initial_ensamble(double *sites, curandState *state, setting *settings){
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	for (int i = 0; i < (int) settings->n; i++) {
		sites[idx * (int) settings->n + i] = rand_x(state + idx, -100, 100);
	}
}

void collectConfigs(config *settings, char **argv){
	settings->runs 			= atof(argv[1]);
	settings->SII 			= atof(argv[2]);
	settings->RS 			= atof(argv[3]);
	settings->xjPerThread 		= atof(argv[4]);
	settings->n 			= atof(argv[5]);
	settings->threadsPerBlock 	= atof(argv[6]);
	settings->MCS 			= atof(argv[7]);
	settings->MCI 			= atof(argv[8]);
	settings->a 			= atof(argv[9]);
	settings->m0 			= atof(argv[10]);
	settings->lambda 		= atof(argv[11]);
	settings->mu_sq 		= atof(argv[12]);
	settings->f_sq 			= atof(argv[13]);

	settings->total_configurations = setting->MCS * setting->runs * setting->MCI / setting->SII;
	settings->array_size = settings->MCS * N * sizeof(double);
	settings->blocksPerGrid = settings->MCS;
	settings->N = settings->xjPerThread * settings->threadsPerBlock;

}


int main(int argc, char** argv) {
	printf("Arguments:");
	for(size_t i = 0; i < argc; i++){
		printf(" %i: %s,", i, argv[i]);
	}
	printf("\n");
	printf("Collecting Settings..\n"); // initialise Settings: data_path, runs, SII, RS, xjPerThread, n, threadsPerBlock, MCS, MCI, a, m0, lambda, mu_sq, f_sq
	config *h_settings, *h_settings;
	gpuErrchk(cudaMallocHost((void **)&h_settings, sizeof(setting)));
	gpuErrchk(cudaMalloc((void **)&h_settings, sizeof(setting)));

	collectConfigs(h_settings, argv);

	gpuErrchk(cudaMemcpy(d_settings, h_settings, sizeof(setting), cudaMemcpyHostToDevice));

	
	//printf("Settings: %i runs, %i SII, %i RS, %i xjPerthread, n=%i, %i threadsPerBlock, %i MCS, %i MCI, a=%lf, m0=%lf, lambda=%lf, mu_sq=%lf, f_sq=%lf\n", runs, SII, RS, xjPerThread, (int) h_settings[4], threadsPerBlock, MCS, MCI, h_settings[8], h_settings[9], h_settings[10], h_settings[11], h_settings[12]);



	printf("Degrees of freedom N: %ld \n", N);

	printf("Launching.\n"); // setup main configs
	FILE *data_file;
	data_file = fopen(argv[1], "w");
	srand(42);


	// setup cuRAND state for randomization
	printf("Setting cuRAND..\n");
	curandState *d_state;
	gpuErrchk(
			cudaMalloc((void **)&d_state, h_config->threadsPerBlock * h_settings->blocksPerGrid * sizeof(curandState)));
	init_rand<<<h_settings->blocksPerGrid, h_conifg->threadsPerBlock>>>(d_state);
	gpuErrchk(cudaDeviceSynchronize());

	// setup kernel configs
	printf("Threads per block: %ld; Blocks per grid: %ld\n", h_settings->threadsPerBlock,
			h_settings->blocksPerGrid);

	// initialise device sites
	double *h_sites, *d_sites;
	gpuErrchk(cudaMallocHost((void **)&h_sites, h_settings->array_size));
	gpuErrchk(cudaMalloc((void **)&d_sites, h_settings->array_size));
	gpuErrchk(cudaMemset(h_sites, 0, h_settings->array_size));
	gpuErrchk(cudaMemset(d_sites, 0, h_settings->array_size));


	clock_t start = clock();
	for (int k = 0; k < runs; k++) {
		gpuErrchk(cudaMemcpy(d_sites, h_sites, h_settings->array_size, cudaMemcpyHostToDevice));
		initial_ensamble<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_sites, d_state, d_settings);

		// wait until the Simulation comes to equilibrium and take first messurement
		Simulate<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_sites, RS, d_state, d_settings);

		size_t iterations = MCI / SII;
		for (size_t i = 0; i < iterations; i++) {
			gpuErrchk(cudaDeviceSynchronize());
			gpuErrchk(cudaMemcpy(h_sites, d_sites, h_settings->array_size, cudaMemcpyDeviceToHost));
			Simulate<<<h_settings->blocksPerGrid, h_settings->threadsPerBlock>>>(d_sites, SII, d_state, d_settings);
			messure(data_file, h_sites, h_settingss);
			printf("(%3.2lf %%)\n",
					(i + 1 + k * iterations) * 100.0 / iterations / runs);
		}
	}
	clock_t end = clock();
	fclose(data_file);
	gpuErrchk(cudaDeviceSynchronize());

	// free memory
	printf("Free memory..\n");
	gpuErrchk(cudaFreeHost(h_sites));
	gpuErrchk(cudaFree(d_sites));

	printf("Done.\n");
	double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
	double benchmark = MCI * MCS * N * runs / time_spend;
	printf("Configurations Messured: %ld; Time Spend on Routine: %0.2lfs; "
			"Benchmark: %0.2e sites/s\n",
			total_configurations, time_spend, benchmark);
	return 0;
}
