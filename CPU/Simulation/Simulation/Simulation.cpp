#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <time.h>

using namespace std;

// Settings
constexpr size_t SITES_COUNT = 1000; // before (4.16) Fig. 4
constexpr size_t MONT_CRLO_ITER =
    50; // N_t = 50 (4.8) --- with exp(-d_action) criteria N_t = 10 is
        // sufficiant from experience
constexpr size_t MONT_CRLO_ALGO = 10;         // Resulting from N_t and N_E
constexpr size_t STATISTICAL_INDEPENDENT_ITER = 5; // before (4.9)
constexpr size_t CONFIGURATIONS =
    MONT_CRLO_ITER * MONT_CRLO_ALGO /
    STATISTICAL_INDEPENDENT_ITER; // N_E = 10^2 before (4.11)

constexpr double MU = 1.0;   // before (4.16) Fig. 4
constexpr double LAMBDA = 0; // before (4.16) Fig. 4

double A = 1;             // before (4.16) Fig. 4
constexpr double M_0 = 1; // before (4.16) Fig. 4
double DELTA = 2 * sqrt(A); // 2*sqrt(a) (4.9)
constexpr int N = 10;       // n_tilde = 10 (4.10)

string FILE_PATH("data_fig5.csv");

// Crystal
double *sites = new double[SITES_COUNT];

static double action() {
  double temp_sum = 0;
  for (size_t j = 0; j < SITES_COUNT - 1; j++)
    temp_sum += M_0 / 2 * pow(sites[j] - sites[j + 1], 2) / A +
                (pow(MU * sites[j], 2) / 2 + LAMBDA * pow(sites[j], 4)) * A;

  return temp_sum;
}

int main() {
  cout << "Generating " << CONFIGURATIONS << " Configurations every "
       << STATISTICAL_INDEPENDENT_ITER << " Iterations..\n";

  // uniform distributor
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(-0.5, 0.5);

  ofstream file;
  file.open(FILE_PATH);

  clock_t start = clock();
  for (size_t k = 0; k < MONT_CRLO_ALGO; k++) {

    // initial ensamble
    for (size_t j = 0; j < SITES_COUNT - 1; j++) {
      sites[j] = dis(gen);
      file << sites[j] << ";";
    }

    // x_0 = x_N
    sites[SITES_COUNT - 1] = sites[0];
    file << sites[SITES_COUNT - 1] << ";";

    // Monte Carlo Iterations
    for (size_t i = 0; i < MONT_CRLO_ITER; i++) {
      cout << "("
           << (k * MONT_CRLO_ITER + i + 1) /
                  (MONT_CRLO_ALGO * MONT_CRLO_ITER * 1.0) * 100
           << "%) " << k + 1 << ". Monte Carlo Algorythm: " << i + 1
           << ". Monte Carlo Iteration: S=" << action() << "\n";
      if (i % STATISTICAL_INDEPENDENT_ITER == 0)
        file << "\n";

      // Markov Chain
      for (size_t j = 0; j < SITES_COUNT; j++) {
        for (size_t _ = 0; _ < N; _++) {
          double position = sites[j];
          double new_position = position + dis(gen);

          double old_action = action();
          sites[j] = new_position;
          double d_action = action() - old_action;
          sites[j] = position;

          double r = dis(gen);

          if (d_action < 0)
            sites[j] = new_position;
          // else if (position - DELTA <= new_position &&
          // position + DELTA >= new_position)
          // sites[j] = new_position;
          else if (pow(M_E, -d_action) > r)
            sites[j] = r;
        }
        if (i % STATISTICAL_INDEPENDENT_ITER == 0)
          file << sites[j] << ";";
      }
    }
  }
  clock_t end = clock();

  file.close();
  delete[] sites;

  double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
  double benchmark = MONT_CRLO_ITER * MONT_CRLO_ALGO * SITES_COUNT / time_spend;
  printf("%6.3fs spend. Benchmark: %6.3f sites/s.\n", time_spend, benchmark);
  return 0;
}
