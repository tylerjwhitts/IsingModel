#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;

// Function for generating a random matrix of size (rows, cols) with up or down spins represented by +1, -1
vector<vector<int>> generateRandomSpinMatrix(size_t rows, size_t cols){

    vector<vector<int>> spins(rows, vector<int>(cols));

    // Seed the random number generator with current time
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator(seed); // Mersenne Twister engine

    // Define a uniform integer distribution between 0 and 1
    uniform_int_distribution<int> distribution(0, 1);

    // Fill spin matrix
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            spins[i][j] = (distribution(generator)==0) ? -1 : 1;
        }
    }
    return spins;
}


double Energy(vector<vector<int>> spins){
    size_t rows = spins.size();
    size_t cols = spins[0].size();
    double total_energy = 0;

    for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols; ++j) {
            int nn_sum = 0;
            vector<vector<int>> nearestNeighbourVecs = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
            for (int vec = 0; vec < nearestNeighbourVecs.size(); ++vec){
                int di = nearestNeighbourVecs[vec][0];
                int dj = nearestNeighbourVecs[vec][1];
                // Find nearest neighbours with periodic boundary conditions
                int ii = ((i + di) + rows ) % rows; 
                int jj = ((j + dj) + cols ) % cols;
                nn_sum += spins[ii][jj];
            }
            total_energy += - (1 / 2) * spins[i][j] * nn_sum;
        }
    }
    return total_energy;
}

int deltaE(vector<vector<int>> spins, vector<int> spin_to_flip){
    int i = spin_to_flip[0];
    int j = spin_to_flip[1];
    size_t rows = spins.size();
    size_t cols = spins[0].size();

    int nn_sum = 0;
    vector<vector<int>> nearestNeighbourVecs = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    for (int vec = 0; vec < nearestNeighbourVecs.size(); ++vec){
        int di = nearestNeighbourVecs[vec][0];
        int dj = nearestNeighbourVecs[vec][1];
        int ii = ((i + di) + rows ) % rows;
        int jj = ((j + dj) + cols ) % cols;
        nn_sum += spins[ii][jj];
    }
    return 2 * spins[i][j] * nn_sum;
}

int main() {

    // int num_sweeps = 10000;
    // int N = 500;
    size_t L = 3;

    vector<vector<int>> spins = generateRandomSpinMatrix(L, L);

    cout << "Random spin matrix is: " << endl;
    for (size_t i = 0; i < L; ++i) {
        for (size_t j = 0; j < L; ++j) {
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }

    int energy = Energy(spins);

    cout << "Total energy of the spin matrix is: " << energy <<endl;

    cout << "Flipping a random spin and calculating the energy difference: " << endl;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(1, L);
    vector<int> spin_to_flip = { dist(gen), dist(gen) };

    cout << " Random spin to flip is at " << spin_to_flip[0] << ", " << spin_to_flip[1] << endl;

    // Calculate deltaE directly
    double dE_direct = deltaE(spins, spin_to_flip);

    // Calculate deltaE manually using Energy function
    int iFlip = spin_to_flip[0];
    int jFlip = spin_to_flip[1];
    spins[iFlip][jFlip] = - spins[iFlip][jFlip];

    cout << "Spin matrix after flip is: " << endl;
    for (size_t i = 0; i < L; ++i) {
        for (size_t j = 0; j < L; ++j) {
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }

    double energy_after_flip = Energy(spins);

    double dE_manual = energy_after_flip - energy;

    cout << "The value of deltaE calculated using the deltaE function is " << dE_direct << endl;
    cout << "The value of deltaE calculated using the Energy function is " << dE_manual << endl;

// for (int sweep=0; sweep<num_sweeps; sweep++){
//     for (int i=0; i<N;i++){

//     }
// }

    return 0;
}