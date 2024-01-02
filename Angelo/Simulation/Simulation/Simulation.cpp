#include <iostream>
#include <random>
#include <fstream>

using namespace std;


// Settings
constexpr size_t SITES_COUNT = 1000; // before (4.16) Fig. 4
constexpr size_t MONT_CRLO_ITER = 50; // N_t = 50 (4.8)
constexpr size_t MONT_CRLO_ALGO = 10; // Resulting from N_t and N_E
constexpr size_t STATISTICAL_INDEPENDENT_ITER = 5; // befor (4.9)
constexpr size_t CONFIGURATIONS = MONT_CRLO_ITER * MONT_CRLO_ALGO / STATISTICAL_INDEPENDENT_ITER; // N_E = 10^2 before (4.11)

constexpr double MU = 1; // before (4.16) Fig. 4
constexpr double LAMBDA = 0; // before (4.16) Fig. 4

double A = 0.1; // before (4.16) Fig. 4
constexpr double M_0 = 1.0; // befor (4.16) Fig. 4
double DELTA = 2.0*sqrt(A); // 2*sqrt(a) (4.9)
constexpr int N = 10; // n_tilde = 10 (4.10)

string FILE_PATH("data.csv");

// Crystal
double* sites = new double[SITES_COUNT];

static double action()
{
    double temp_sum = 0;
    for (size_t j = 0; j < SITES_COUNT-1; j++)
        temp_sum += M_0 / 2 * pow(sites[j] - sites[j+1], 2) / A + (pow(MU*sites[j], 2) / 2 + LAMBDA * pow(sites[j], 4)) * A;

    return temp_sum;
}

int main()
{
    cout << "Generating " << CONFIGURATIONS << " Configurations every " << STATISTICAL_INDEPENDENT_ITER << " Iterations..\n";

    // uniform distributor
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    ofstream file;
    file.open(FILE_PATH);

    for (size_t k = 0; k < MONT_CRLO_ALGO; k++) {

        // initial ensamble
        for (size_t j = 0; j < SITES_COUNT-1; j++) {
            sites[j] = dis(gen);
            file << sites[j] << ";";
        }

        // x_0 = x_N
        sites[SITES_COUNT - 1] = sites[0];
        file << sites[SITES_COUNT - 1] << ";";

        // Monte Carlo Iterations
        for (size_t i = 0; i < MONT_CRLO_ITER; i++) {
            cout << "(" << (i + 1) * (k + 1) / (MONT_CRLO_ALGO * MONT_CRLO_ITER * 1.0) * 100 << "%) " << k << ". Monte Carlo Algorythm: " << i << ". Monte Carlo Iteration\n";
            if (i%STATISTICAL_INDEPENDENT_ITER==0)
                file << "\n";

            // Markov Chain
            for (size_t j = 0; j < SITES_COUNT; j++) {
                for (size_t _ = 0; _ < N; _++) {
                    double position = sites[j];
                    double new_position = dis(gen);

                    double old_action = action();
                    sites[j] = new_position;
                    double d_action = action() - old_action;
                    sites[j] = position;

                    if (d_action < 0 || (position - DELTA <= new_position && position + DELTA >= new_position))
                        sites[j] = new_position;
                }
                if (i%STATISTICAL_INDEPENDENT_ITER==0)
                    file << sites[j] << ";";
            }
        }

    }
    file.close();
    delete[] sites;

    return 0;
}