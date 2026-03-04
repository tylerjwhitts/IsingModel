#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <pybind11/pybind11.h>

using namespace std;

namespace py = pybind11;

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


float Energy(vector<vector<int>> spins){
    size_t rows = spins.size();
    size_t cols = spins[0].size();
    float total_energy = 0;

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
            total_energy += - 0.5 * spins[i][j] * nn_sum;
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

void printSpinMatrix(vector<vector<int>> spins){
    size_t rows = spins.size();
    size_t cols = spins[0].size();
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }
}

void runEnergyCalculationTests(size_t latticeLength) {

    vector<vector<int>> spins = generateRandomSpinMatrix(latticeLength, latticeLength);

    cout << "Random spin matrix is: " << endl;
    printSpinMatrix(spins);

    int energy = Energy(spins);

    cout << "Total energy of the spin matrix is: " << energy << "\n";

    cout << "Flipping a random spin and calculating the energy difference: " << endl;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, latticeLength-1);
    vector<int> spin_to_flip = { dist(gen), dist(gen) };

    cout << "Random spin to flip is at (" << spin_to_flip[0] << ", " << spin_to_flip[1] << ")" << endl;

    // Calculate deltaE directly
    int dE_direct = deltaE(spins, spin_to_flip);

    // Calculate deltaE manually using Energy function
    int iFlip = spin_to_flip[0];
    int jFlip = spin_to_flip[1];
    spins[iFlip][jFlip] = - spins[iFlip][jFlip];

    cout << "Spin matrix after flip is: " << endl;
    printSpinMatrix(spins);

    int energy_after_flip = Energy(spins);

    int dE_manual = energy_after_flip - energy;

    cout << "The value of deltaE calculated using the deltaE function is " << dE_direct << endl;
    cout << "The value of deltaE calculated using the Energy function is " << dE_manual << endl;
}

void runMarkovChain(size_t L, double temperature, int numSteps, int numSweeps, int transientSweeps=20, bool outputProgress=false){
    
    vector<vector<int>> spins = generateRandomSpinMatrix(L, L);
    size_t rows = spins.size();
    size_t cols = spins[0].size();

    // Random number generator for spin coordinates to flip
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> spinCoordDist(0, L-1);

    // Random number generator for spin for acceptance criterion
    uniform_real_distribution<double> probDist(0.0, 1.0);

    cout << " Initial spin configuration: " << endl;
    printSpinMatrix(spins);
    double initialEnergy = Energy(spins) / (L * L);
    cout << "The energy of the initial configuration is: " << initialEnergy << endl;

    for (int sweep=0; sweep<numSweeps; sweep++){
        for (int step=0; step<numSteps; step++){
            vector<int> spin_to_flip = { spinCoordDist(gen), spinCoordDist(gen) };
            int iFlip = spin_to_flip[0];
            int jFlip = spin_to_flip[1];
            double dE = deltaE(spins, spin_to_flip);
            double exponent = - dE / temperature;
            double acceptanceProb = exp(exponent);
            if (dE < 0) {
                spins[iFlip][jFlip] = - spins[iFlip][jFlip];
            }
            else if (probDist(gen) < acceptanceProb) {
                spins[iFlip][jFlip] = - spins[iFlip][jFlip];
            }
        }

        if ( outputProgress == true ){
            cout << "Sweep " << sweep+1 << " of " << numSweeps << " complete." << endl;
        }
    }

    cout << "Final spin configuration after " << numSweeps << " sweeps: " << endl;
    printSpinMatrix(spins);

    double finalEnergy = Energy(spins) / (L * L);
    cout << "The normalised energy of final configuration is: " << finalEnergy << endl;

}

PYBIND11_MODULE(IsingMarkovChainMonteCarlo, handle) {
    handle.doc() = R"(The module runs a markov chain for the monte carlo simulation. 
                      On a specific ising spin lattice at a specified temperature and number of sweeps.
                      Observables energy and magnetisation are recorded at each sweep.)";

    handle.def("runIsingMarkovChainCPP", &runMarkovChain, "Runs the MarkovChain on an Ising Model spin lattice.",
                py::arg("L"), py::arg("T"), py::arg("numSteps"), py::arg("numSweeps"), py::arg("transientSweeps")=20, py::arg("outputProgress")=true); 
}


int main() {

    int numSweeps = 1000;
    int numSteps = 100;
    int transient_sweeps = 20;
    
    size_t L;
    double temperature;

    cout << "Input desired square lattice length: ";
    cin >> L;
    cout << "\n";

    cout << "Input desired temperature in units of T_0: ";
    cin >> temperature;
    cout << "\n";

    bool outputProgress = false;

    // Initial testing
    // runEnergyCalculationTests(L);

    runMarkovChain(L, temperature, numSteps, numSweeps);

    return 0;
}