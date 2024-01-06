#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <stdio.h>
#include <time.h>

using namespace std;

// Settings
constexpr size_t SITES_COUNT = 800; // (N) before (4.16) Fig. 4
constexpr size_t MONT_CRLO_ITER =
    50; // N_t = 50 (4.8) --- with exp(-d_action) criteria N_t = 10 is
        // sufficiant from experience
constexpr size_t MONT_CRLO_ALGO = 100;              // Resulting from N_t and N_E
constexpr size_t STATISTICAL_INDEPENDENT_ITER = 5; // before (4.9)
int RECORDING_START = 0;
constexpr size_t CONFIGURATIONS =
    MONT_CRLO_ITER * MONT_CRLO_ALGO /
    STATISTICAL_INDEPENDENT_ITER; // N_E = 10^2 before (4.11)

double MU_SQ = 1.0;          // before (4.16) Fig. 4
constexpr double LAMBDA = 0; // before (4.16) Fig. 4
double ff = 2;

double A = 1.0;             // before (4.16) Fig. 4
constexpr double M_0 = 1.0; // before (4.16) Fig. 4
double DELTA = 2 * sqrt(A); // 2*sqrt(a) (4.9)
constexpr int N = 10;       // n_tilde = 10 (4.10)

string FILE_PATH("data_test.csv");
constexpr double initial_ensamble_radius = 1;

int actualMessurements = 0;

// Crystal
double *sites = new double[SITES_COUNT];

static double Potential1(double *xj) {
  return (MU_SQ * pow(*xj, 2) / 2 + LAMBDA * pow(*xj, 4));
}

static double Potential2(double *xj) {
  return LAMBDA * pow(pow(*xj, 2) - ff, 2);
}

static double getTotalAction() {
  double temp_sum = 0;
  for (size_t j = 0; j < SITES_COUNT - 1; j++)
    temp_sum += M_0 / 2 * pow(sites[j] - sites[j + 1], 2) / A +
                A * Potential1(&sites[j]);
  return temp_sum;
}

static double getAction(double *xj) {
  return M_0 / 2 * pow(*xj - *(xj + 1), 2) / A + A * Potential1(xj);
}

static double getDAction(double *xj, double *new_xj) {
  double old_action = getAction(xj) + getAction(xj - 1);
  double temp = *xj;
  *xj = *new_xj;
  double new_action = getAction(xj) + getAction(xj - 1);
  *xj = temp;
  return new_action - old_action;
}

void messure(ofstream *file) {
  for (int i = 0; i < SITES_COUNT - 1; i++) {
    *file << sites[i] << ";";
  }
  *file << sites[SITES_COUNT - 1] << "\n";
}

static double rand_x(double lowerBound, double upperBound) {
  return double(rand()) / double(RAND_MAX) * (upperBound - lowerBound) +
         lowerBound;
}

int main() {
  printf("Generating %lu configurations every %lu Iterations..\n",
         CONFIGURATIONS, STATISTICAL_INDEPENDENT_ITER);

  srand(42);

  ofstream file;
  file.open(FILE_PATH);
  for (int i = 0; i < SITES_COUNT - 1; i++)
    file << i << ";";
  file << SITES_COUNT - 1 << "\n";

  clock_t start = clock();
  for (size_t k = 0; k < MONT_CRLO_ALGO; k++) {

    // initial ensamble
    for (size_t j = 0; j < SITES_COUNT - 1; j++) {
      sites[j] = rand_x(initial_ensamble_radius, -initial_ensamble_radius);
    }

    // x_0 = x_N
    sites[SITES_COUNT - 1] = sites[0];

    // messure(&file);

    // Monte Carlo Iterations
    for (size_t i = 0; i < MONT_CRLO_ITER; i++) {
      if (false) {
        printf("(%3.2f %%) %3d/%3d (MCS/MCI): S=%0.2f\n",
               double((k * MONT_CRLO_ITER + i + 1) /
                      (MONT_CRLO_ALGO * MONT_CRLO_ITER * 1.0) * 100),
               int(k + 1), int(i + 1), getTotalAction());
      } else {

        printf("(%3.2f %%)\n",
               double((k * MONT_CRLO_ITER + i + 1) /
                      (MONT_CRLO_ALGO * MONT_CRLO_ITER * 1.0) * 100));
      }

      // Markov Chain
      for (size_t j = 1; j < SITES_COUNT - 1; j++) {
        for (size_t _ = 0; _ < N; _++) {
          double new_xj = rand_x(sites[j] - DELTA, sites[j] + DELTA);
          double d_action = getDAction(&sites[j], &new_xj);

          if (d_action < 0)
            sites[j] = new_xj;
          else if (pow(M_E, -d_action) > rand_x(0, 1))
            sites[j] = new_xj;
        }
      }
      if (i >= RECORDING_START && (i + 1) % STATISTICAL_INDEPENDENT_ITER == 0)
        messure(&file);
    }
  }
  clock_t end = clock();

  file.close();
  delete[] sites;

  double time_spend = (double)(end - start) / CLOCKS_PER_SEC;
  double benchmark = MONT_CRLO_ITER * MONT_CRLO_ALGO * SITES_COUNT / time_spend;
  printf("%6.3fs spend. Benchmark: %e sites/s.\n", time_spend, benchmark);
  printf("%d", actualMessurements);
  return 0;
}
