#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <list>
#include <tracy/Tracy.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace std;

// Function for generating a random matrix of size (rows, cols) with up or down spins represented by +1, -1
vector<vector<int>> generateRandomSpinMatrix(size_t rows, size_t cols){
    ZoneScoped;
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


float Energy(const vector<vector<int>> &spins){
    ZoneScoped;
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

int deltaE(const vector<vector<int>> &spins, const vector<int> &spin_to_flip){
    ZoneScoped;
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

double getSpinSum(const vector<vector<int>> &spins){
    size_t rows = spins.size();
    size_t cols = (rows > 0) ? spins[0].size() : 0;
    double spinSum = 0.0;
    for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols; ++j) {
            spinSum += spins[i][j];
        }
    }
    return spinSum;
}

double getMagSq(const vector<vector<int>> &spins){
    ZoneScoped;
    size_t rows = spins.size();
    size_t cols = (rows > 0) ? spins[0].size() : 0;
    size_t total_size = rows * cols;
    double spin_sum = 0.0;
    for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols; ++j) {
            spin_sum += spins[i][j];
        }
    }
    double mag = spin_sum / total_size;
    return mag * mag;
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

tuple<vector<double>, vector<double>> runMarkovChain(size_t L, double temperature, int numSteps, int numSweeps, int transientSweeps=20, bool completeMeasure=false, bool outputProgress=false){
    
    ZoneScopedN("MarkovChainAlgorithm")

    vector<vector<int>> spins = generateRandomSpinMatrix(L, L);
    size_t rows = spins.size();
    size_t cols = spins[0].size();

    // Random number generator for spin coordinates to flip
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> spinCoordDist(0, L-1);

    // Random number generator for spin for acceptance criterion
    uniform_real_distribution<double> probDist(0.0, 1.0);

    // Initial energy of the system and total spin sum
    double isingEnergy = Energy(spins);
    double spinSum = getSpinSum(spins);

    // Initial spin sum


    if (outputProgress == true) {
        cout << " Initial spin configuration: " << endl;
        printSpinMatrix(spins);
        cout << "The energy of the initial configuration is: " << isingEnergy << endl;
    }

    // Initialise lists for observable measurements
    vector<double> energyMeasurements = {};
    vector<double> magSqMeasurements = {};

    for (int sweep=0; sweep<numSweeps+transientSweeps; sweep++){
        for (int step=0; step<numSteps; step++){
            vector<int> spin_to_flip = { spinCoordDist(gen), spinCoordDist(gen) };
            int iFlip = spin_to_flip[0];
            int jFlip = spin_to_flip[1];
            double dE = deltaE(spins, spin_to_flip);
            double acceptanceProb = exp(- dE / temperature);
            if (dE < 0) {
                spins[iFlip][jFlip] = - spins[iFlip][jFlip];
                isingEnergy += dE;
                spinSum += 2*spins[iFlip][jFlip]; // Add spin twice to remove old spin
            }
            else if (probDist(gen) < acceptanceProb) {
                spins[iFlip][jFlip] = - spins[iFlip][jFlip];
                isingEnergy += dE;
                spinSum += 2*spins[iFlip][jFlip];
            }
            FrameMark;
        }

        // Measure
        if (sweep >= transientSweeps){
            if (completeMeasure==true){
            energyMeasurements.push_back( Energy(spins) / (L*L) );
            magSqMeasurements.push_back(getMagSq(spins));
            }
            else {
            // Switch to storing the energy and mag of system 
            // and altering it by deltaE or -2S each time a spin flip is accepted
            // Complexity goes from O(N) to O(1)
            energyMeasurements.push_back( isingEnergy / (L * L));
            double mag = spinSum / (L * L);
            magSqMeasurements.push_back(mag * mag);
            }
        }

        if ( outputProgress == true ){
            cout << "Sweep " << sweep+1 << " of " << numSweeps << " complete." << endl;
        }
        FrameMark;
    }

    if (outputProgress == true) {
        cout << "Final spin configuration after " << numSweeps << " sweeps: " << endl;
        printSpinMatrix(spins);

        double finalEnergy = Energy(spins) / (L * L);
        cout << "The normalised energy of final configuration is: " << finalEnergy << endl;
    }

    return make_tuple(energyMeasurements, magSqMeasurements);
}

PYBIND11_MODULE(IsingMarkovChainMonteCarlo, handle) {
    handle.doc() = R"(The module runs a markov chain for the monte carlo simulation. 
                      On a specific ising spin lattice at a specified temperature and number of sweeps.
                      Observables energy and magnetisation are recorded at each sweep.)";

    handle.def("runIsingMarkovChainCPP", &runMarkovChain, "Runs the MarkovChain on an Ising Model spin lattice.",
                py::arg("L"), py::arg("T"), py::arg("numSteps"), py::arg("numSweeps"), py::arg("transientSweeps")=20, py::arg("completeMeasure")=true, py::arg("outputProgress")=true); 
}


double expectationValue(vector<double> &observableMeasurements){
    if (observableMeasurements.empty()){
        return 0.0;
    }

    double sum = accumulate(observableMeasurements.begin(), observableMeasurements.end(), 0);

    return sum / observableMeasurements.size();
}


int main() {

    int numSweeps = 1000;
    int numSteps = 100;
    int transientSweeps = 20;
    
    size_t L;
    double temperature;

    cout << "Input desired square lattice length: ";
    cin >> L;
    cout << "\n";

    cout << "Input desired temperature in units of T_0: ";
    cin >> temperature;
    cout << "\n";

    // L = 3;
    // temperature = 2.0;

    // Initial testing
    // runEnergyCalculationTests(L);

    // Testing the runMarkovChainFunction
    
    // Without completely measuring observables each time
    vector<double> energyMeasurements;
    vector<double> magSqMeasurements; 

    auto start = chrono::steady_clock::now();
    tie(energyMeasurements, magSqMeasurements) = runMarkovChain(L, temperature, numSteps, numSweeps);
    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start);

    double energyExpectation = expectationValue(energyMeasurements);
    double magSqExpectation = expectationValue(magSqMeasurements);
    

    cout << "The Markov Chain algorithm took " << duration.count() << " seconds without complete measurements each time." << endl;
    cout << "The expectation values for the energy and magnetisation squared are, respectively: " << energyExpectation << " and " << magSqExpectation << endl;

    // Completely measuring observables each time
    bool completeMeasure = true;
    vector<double> energyMeasurementsComplete;
    vector<double> magSqMeasurementsComplete; 
    auto startComplete = chrono::steady_clock::now();
    tie(energyMeasurementsComplete, magSqMeasurementsComplete) = runMarkovChain(L, temperature, numSteps, numSweeps, transientSweeps, completeMeasure);
    auto endComplete = chrono::steady_clock::now();

    auto durationComplete = chrono::duration_cast<chrono::seconds>(endComplete - startComplete);

    double energyExpectationComplete = expectationValue(energyMeasurementsComplete);
    double magSqExpectationComplete = expectationValue(magSqMeasurementsComplete);
    

    cout << "The Markov Chain algorithm took " << durationComplete.count() << " seconds when performing complete measurements each time." << endl;
    cout << "The expectation values for the energy and magnetisation squared are, respectively: " << energyExpectationComplete << " and " << magSqExpectationComplete << endl;

    return 0;
}