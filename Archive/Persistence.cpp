//? Changed jsonData.size()/2.0 into jsonData.size()-1.0
//? Calibrated the gap values to actual time values using timestep and dataCompression
#include <iostream>
#include <fstream>
#include <vector>
#include "json.hpp"
#include <cmath>

using namespace std;

void persistence(int numParticles){
    int dataCompression=500;
    int numSteps = 100000;
    double timestep=0.001;
    //load the tensor data
    std::ifstream file("particle_positions.json");
    if (!file.is_open()) {
        cerr << "Unable to open input file." << endl;
    }

    nlohmann::json jsonData;
    file >> jsonData;
    cout << "Number of gaps: " << jsonData.size() << endl;
    file.close();

    //create an empty list to store the persistence list of each partilce at each gap to take average later
    vector<vector<double>> persistencePerParticle(numParticles,vector<double>(floor(jsonData.size()/2.0),0.0));
    //create a for loop running through each particle
    for (int i = 0; i < numParticles; ++i) {
        //create an empty list to store the persistence values of each gap
        vector<double> persistencePerGap(floor(jsonData.size()-1.0), 0.0);
        //create a for loop running through half the size of the number of timesteps; here j will be caled "gap"
        for (int j = 0; j < floor(jsonData.size()-1.0); ++j) {
            for (int k = 0; k < (jsonData.size()-j-1); ++k) {
                // for each k calculate (cos(jsonData[k][i][4])*cos(jsonData[k+j+1][i][4]))/(cos(jsonData[k][i][4]*cos(jsonData[k][i][4]) and add it to the list
                // cos function is being called with an argument of type nlohmann::json::value_type which is a double
                persistencePerGap[j] += (cos(static_cast<double>(jsonData[k][i][4])) * cos(static_cast<double>(jsonData[k+j+1][i][4])))  ;
            }
            //take the average of the list and store it in the list
            persistencePerGap[j] /= (jsonData.size()-j-1);
        }
        //Add the list to the list of lists
        persistencePerParticle[i]=persistencePerGap;
    }

    // create a list which stores the average persistence of each particle at a given gap
    vector<double> averagePersistencePerGap(floor(jsonData.size()-1.0), 0.0);
    //create a for loop running through half the size of the number of timesteps; here j will be caled "gap"
    for (int j = 0; j < floor(jsonData.size()-1.0); ++j) {
        //create a for loop running through each particle
        for (int i = 0; i < numParticles; ++i) {
            //add the persistence of each particle at a given gap to the list
            averagePersistencePerGap[j] += persistencePerParticle[i][j];
        }
        //take the average of the list and store it in the list
        averagePersistencePerGap[j] /= numParticles;
    }

    // Plot averagePersistencePerGap vs gap
    ofstream outputFile("persistence.dat");
    if (!outputFile.is_open()) {
        cerr << "Unable to open output file." << endl;
    }
    for (int i = 0; i < floor(jsonData.size()-1.0); ++i) {
        outputFile << i*dataCompression*timestep << " " << averagePersistencePerGap[i] << endl; 
    }
    outputFile.close();

    //Plot persistence.dat using gnu plot
    system("gnuplot -p -e \"set xlabel 'Gap'; set ylabel 'Persistence'; plot 'persistence.dat' using 1:2 with linespoints\"");

};

int main() {
    persistence(100);
    return 0;
}
