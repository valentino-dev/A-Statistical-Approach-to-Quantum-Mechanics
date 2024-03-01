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

typedef struct{
	int runs, SII, RS, xjPerThread, n, threadsPerBlock, MCS, MCI, a, m0, lambda, mu_sq, f_sq, delta, N, total_configurations, blocksPerGrid, array_size;
} setting;

double rand_x(double min, double max) {
	// return 0.5;
	return rand()/RAND_MAX * (max - min) + min;
}

double potential_1(double *x, setting *settings) {
	return settings->mu_sq * pow(*x, 2) * 0.5 + settings->lambda * pow(*x, 4);
}

double potential_2(double *x, setting *settings) {
	return settings->lambda * pow(*x * *x - settings->f_sq, 2);
}

double calc_S_of_xj(double x, double *ptr, double *previous_site, double *following_site, setting *settings) {
	return settings->m0 * 0.5 * (pow(*following_site - x, 2) + pow(x - *previous_site, 2)) / settings->a +
		settings->a * (potential_1(&x, settings) + potential_1(ptr - 1, settings));
}

double calc_dS(double *xptr, double *previous_site, double *following_site, double *newx_ptr, setting *settings) {
	return calc_S_of_xj(*newx_ptr, xptr, previous_site, following_site, settings) - 
		calc_S_of_xj(*xptr, xptr, previous_site, following_site, settings);
}

void step(double *sites, int site, int previous_site, int following_site, setting *settings) {
	double *xptr = sites+site;
	double new_x = rand_x(*xptr - settings->delta, *xptr + settings->delta);
	double dS = calc_dS(xptr, sites+previous_site, sites+following_site, &new_x, settings);
	if (dS < 0 || pow(M_E, -dS) > rand_x(0.0, 1.0))
		*xptr = new_x;
}

void Print_Action(double *sites, setting *settings){
	double S = 0;
	for (int i = 0; i < settings->N; i++)
		S += settings->m0 * 0.5 * (pow(sites[(i+1)%settings->N] - sites[i], 2)) / settings->a + settings->a * potential_1(&sites[i], settings);

	printf("S: %lf\n", S);
}

// cuda c kernel
void Simulate(double *sites, int iterations, setting *settings) {
	for (int _ = 0; _ < iterations; _++) {
		for (int i = 0; i < settings->N)
				int site = i;
				int previous_site = (i - 1) % settings->N;
				int following_site = (i + 1) % settings->N;
				for (int __ = 0; __ < settings->n; __++) {
					step(sites, site, previous_site, following_site, settings);
				}
		}
}

void print_sites(double *sites) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
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

void initial_ensamble(double *sites, setting *settings){
	for (int i = 0; i < (int) settings->n; i++) {
		sites[(int) settings->n + i] = rand_x(state + idx, -100, 100);
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
	h_settings = malloc((void **)&h_settings, sizeof(setting));
	collectSettings(h_settings, argv);
	printf("Settings: %i runs, %i SII, %i RS, %i xjPerthread, n=%i, %i threadsPerBlock, %i MCS, %i MCI, a=%lf, m0=%lf, lambda=%lf, mu_sq=%lf, f_sq=%lf\n", h_settings->runs, h_settings->SII, h_settings->RS, h_settings->xjPerThread, h_settings->n, h_settings->threadsPerBlock, h_settings->MCS, h_settings->MCI, h_settings->a, h_settings->m0, h_settings->lambda, h_settings->mu_sq, h_settings->f_sq);
	printf("Degrees of freedom N: %ld \n", h_settings->N);
	printf("Threads per block: %d; Blocks per grid: %d\n", h_settings->threadsPerBlock, h_settings->blocksPerGrid);


	// Setup main settings
	printf("Launching.\n"); 
	FILE *data_file;
	data_file = fopen(argv[1], "w");

	// Initialise random
	srand(42);

	// Initialise device sites
	double *h_sites, *d_sites;
	h_sites = malloc((void **)&h_sites, h_settings->array_size);


	// Starting runs
	clock_t start = clock();
	for (int k = 0; k < h_settings->runs; k++) {
		initial_ensamble(h_sites, h_settings);


		// Wait until the Simulation comes to equilibrium
		Simulate(h_sites, h_settings->RS, h_settings);


		// Simulating and messuring
		size_t iterations = h_settings->MCI / h_settings->SII;
		for (size_t i = 0; i < iterations; i++) {
			Simulate(h_sites, iterations, h_settings);
			messure(data_file, h_sites, h_settings);
			printf("(%3.2lf %%)\n", (i + 1 + k * iterations) * 100.0 / iterations / h_settings->runs);
		}
	}
	clock_t end = clock();
	fclose(data_file);


	// Free memory
	printf("Free memory..\n");
	free(h_sites);


	// End Results
	printf("Done.\n");
	double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
	double benchmark = h_settings->MCI * h_settings->MCS * h_settings->N * h_settings->runs / time_spend;
	printf("Configurations Messured: %ld; Time Spend on Routine: %0.2lfs; "
			"Benchmark: %0.2e sites/s\n",
			h_settings->total_configurations, time_spend, benchmark);
	return 0;
}
