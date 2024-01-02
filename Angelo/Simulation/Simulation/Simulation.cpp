#include <iostream>
#include <random>
#include <fstream>

using namespace std;


// Settings
constexpr size_t SITES_COUNT = 1000;
constexpr size_t MONT_CRLO_ITER = 5;

constexpr double MU = 1;
constexpr double LAMBDA = 0;

double A = 0.1;
constexpr double M_0 = 1.0;
double DELTA = 2.0*sqrt(A);
constexpr int N = 10;

string FILE_PATH("data.csv");

// Crystal
double* sites = new double[SITES_COUNT];

static double action()
{
    double temp_sum = 0;
    for (size_t j = 0; j < SITES_COUNT-1; j++)
        temp_sum += M_0 / 2 * pow(sites[j] - sites[j+1], 2) / A + (pow(MU, 2) / 2 * pow(sites[j], 2) + LAMBDA * pow(sites[j], 4)) * A;

    return temp_sum;
}

int main()
{
    // uniform distributor
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    ofstream file;
    file.open(FILE_PATH);

    // initial ensamble
    cout << "Initial ensamble: \n";
    for (size_t j = 0; j < SITES_COUNT-1; j++) {
        sites[j] = dis(gen);
        file << sites[j] << ";";
    }

    // x_0 = x_N
    sites[SITES_COUNT - 1] = sites[0];
    file << sites[SITES_COUNT - 1] << ";";

    // Monte Carlo Iterations
    for (size_t i = 0; i < MONT_CRLO_ITER; i++) {
        cout << "Monte Carlo Iteration: " << i << "\n";
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
            file << sites[j] << ";";
        }
    }

    file.close();
    delete[] sites;

    return 0;
}