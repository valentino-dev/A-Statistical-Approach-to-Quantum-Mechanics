#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>
#include <time.h>

using namespace std;

// Settings
constexpr size_t SITES_COUNT = 1000;               // before (4.16) Fig. 4
constexpr size_t MONT_CRLO_ITER = 50;              // N_t = 50 (4.8)
constexpr size_t SIMULATIONS_COUNT = 10;           // Resulting from N_t and N_E
constexpr size_t STATISTICAL_INDEPENDENT_ITER = 5; // befor (4.9)
constexpr size_t CONFIGURATIONS =
    MONT_CRLO_ITER * SIMULATIONS_COUNT /
    STATISTICAL_INDEPENDENT_ITER; // N_E = 10^2 before (4.11)

constexpr double MU = 1;     // before (4.16) Fig. 4
constexpr double LAMBDA = 0; // before (4.16) Fig. 4

double A = 0.1;               // before (4.16) Fig. 4
constexpr double M_0 = 1.0;   // befor (4.16) Fig. 4
double DELTA = 2.0 * sqrt(A); // 2*sqrt(a) (4.9)
constexpr int N = 10;         // n_tilde = 10 (4.10)

string FILE_PATH("data.csv");

// Crystal
double *sites = new double[SITES_COUNT];

__global__ void MonteCarloSimulation() {

  // initial ensamble
  for (size_t j = 0; j < SITES_COUNT - 1; j++) {
    sites[j] = dis(gen);
    // file << sites[j] << ";";
  }

  // x_0 = x_N
  sites[SITES_COUNT - 1] = sites[0];
  // file << sites[SITES_COUNT - 1] << ";";

  // Monte Carlo Iterations
  for (size_t i = 0; i < MONT_CRLO_ITER; i++) {
    // float percentage = (k * MONT_CRLO_ITER + i + 1) /
    //(SIMULATIONS_COUNT * MONT_CRLO_ITER * 1.0) * 100;
    // printf("(%.2f%%) %i. Monte Carlo Algorythm: %i. Monte Carlo Iteration\n",
    // percentage, k, i); if (i%STATISTICAL_INDEPENDENT_ITER==0) file << "\n";

    // Markov Chain
    for (size_t j = 0; j < SITES_COUNT; j++) {
      for (size_t _ = 0; _ < N; _++) {
        double position = sites[j];

        double new_position = dis(gen);

        double old_action = 0;
        for (size_t j = 0; j < SITES_COUNT - 1; j++)
          old_action +=
              M_0 / 2 * pow(sites[j] - sites[j + 1], 2) / A +
              (pow(MU * sites[j], 2) / 2 + LAMBDA * pow(sites[j], 4)) * A;

        sites[j] = new_position;
        double new_action = 0;
        for (size_t j = 0; j < SITES_COUNT - 1; j++)
          new_action +=
              M_0 / 2 * pow(sites[j] - sites[j + 1], 2) / A +
              (pow(MU * sites[j], 2) / 2 + LAMBDA * pow(sites[j], 4)) * A;
        double d_action = new_action - old_action;

        sites[j] = position;

        if (d_action < 0 || (position - DELTA <= new_position &&
                             position + DELTA >= new_position))
          sites[j] = new_position;
      }
      // if (i%STATISTICAL_INDEPENDENT_ITER==0) file << sites[j] << ";";
    }
  }
}

int main() {
  printf("Generating %i every %i Iterations..\n", CONFIGURATIONS,
         STATISTICAL_INDEPENDENT_ITER);
  double *h_Crystals[SIMULATIONS_COUNT][SITES_COUNT];
  cudaMalloc(&h_Crystals, sizeof(h_Crystals));
  cudaMallocPitch((void **)&array, &pitch, SITES_COUNT * sizeof(double),
                  SIMULATIONS_COUNT);

  ofstream file;
  file.open(FILE_PATH);

  clock_t start = clock();
  MonteCarloSimulation<<<SIMULATIONS_COUNT, 1>>>();

  clock_t end = clock();
  file.close();
  delete[] sites;
  double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
  printf("%6.3fs spend.", time_spend);

  return 0;
}
