#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib> // Include <cstdlib> for std::rand()
#include <numeric> // Include <numeric> for std::accumulate()
#include <fstream> // Include <fstream> for file operations
#include <ctime>
#include <random>

// Particle structure
struct Particle {
    double x, y, vx, vy, ax, ay;
};

// Function to calculate the inter-particle forces
void calculateForces(std::vector<Particle>& particles, double L) {
    // Reset forces
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
    }

    // Calculate forces
    for (int i = 0; i < particles.size() - 1; ++i) {
        for (int j = i + 1; j < particles.size(); ++j) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double r = std::sqrt(dx * dx + dy * dy);

            if (r < 1.0) { // Hard sphere interaction
                double force = 1.0 / r;
                double fx = force * dx/r;
                double fy = force * dy/r;
                particles[i].ax += fx;
                particles[i].ay += fy;
                particles[j].ax -= fx;
                particles[j].ay -= fy;
            }
        }
    }
}

// Function to update the positions and velocities using the Velocity Verlet method
void updatePositionsAndVelocities(std::vector<Particle>& particles, double dt, double L) {
    for (auto& particle : particles) {
        particle.vx += particle.ax * dt;
        particle.vy += particle.ay * dt;
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;

        // Apply periodic boundary conditions
        particle.x -= std::floor(particle.x / L) * L;
        particle.y -= std::floor(particle.y / L) * L;
    }
}

// Function to calculate the total energy of the system
double calculateTotalEnergy(const std::vector<Particle>& particles) {
    double totalEnergy = 0.0;
    for (const auto& particle : particles) {
        double vx = particle.vx;
        double vy = particle.vy;
        totalEnergy += 0.5 * (vx * vx + vy * vy); // Kinetic energy
    }
    return totalEnergy;
}

// Function to calculate the RMS energy fluctuations for various timesteps
void calculateRMSEnergyFluctuations( int NRepeats, int numRuns, const std::vector<double>& timesteps) {
    std::vector<double> rmsEnergyFluctuations(numRuns, 0.0);


    for (int run = 0; run < numRuns; ++run) {

        // Simulation parameters
        int numParticles = 100  ;     // Number of particles
        double boxSize = 100.0;      // Size of the simulation box
        double timestep = timesteps[run];
        int numSteps = std::floor(1000/timestep);        // Number of simulation steps
        
        std::vector<double> energies(NRepeats, 0.0);

        for (int rerun = 0; rerun < NRepeats ; ++rerun) {

            // Initialize particles
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> dis(0.0, 1.0);
            std::vector<Particle> particles(numParticles);
            for (auto& particle : particles) {
                double randomValue1 = dis(gen);
                double randomValue2 = dis(gen);
                particle.x = boxSize * randomValue1 ;
                particle.y = boxSize * randomValue2 ;
                particle.vx = 0.0;
                particle.vy = 0.0;
                particle.ax = 0.0;
                particle.ay = 0.0;
            }


            // Main simulation loop
            for (int step = 0; step < numSteps; ++step) {
                // Calculate forces
                calculateForces(particles, boxSize);

                // Update positions and velocities
                updatePositionsAndVelocities(particles, timestep, boxSize);

            }

            // Calculate and store the energy
            energies[rerun] = calculateTotalEnergy(particles);

            // Calculate RMS energy fluctuation
            // Write energy data to a file for Gnuplot
            // std::ofstream outputFile("energy_data.dat");
            // if (outputFile.is_open()) {
            //     for (int i = 0; i < numSteps; ++i) {
            //         outputFile << i << " " << energies[i] << std::endl;
            //     }
            //     outputFile.close();
            // } else {
            //     std::cerr << "Unable to open output file." << std::endl;
            // }
        
        }
        double averageEnergy = std::accumulate(energies.begin(), energies.end(), 0.0) / NRepeats;
        double sumSquaredDeviations = 0.0;
        for (const auto& energy : energies) {
            double deviation = energy - averageEnergy;
            sumSquaredDeviations += deviation * deviation;
        }
        rmsEnergyFluctuations[run] = std::sqrt(sumSquaredDeviations / NRepeats);
    }
    // Write RMS fluctuation energy data to a file for Gnuplot
    std::ofstream outputFile("rms_energy_fluctuations.dat");
    if (outputFile.is_open()) {
        for (int i= 0; i < numRuns; ++i) {
            outputFile << timesteps[i] << " " << rmsEnergyFluctuations[i] << std::endl;
        }
        outputFile.close();
    } else {
        std::cerr << "Unable to open output file." << std::endl;
    }

    // Gnuplot snippet to plot RMS energy fluctuations
    std::ofstream gnuplotScript("gnuplot_script.plt");
    if (gnuplotScript.is_open()) {
        gnuplotScript << "set xlabel 'Time Step'" << std::endl;
        gnuplotScript << "set ylabel 'RMS Energy Fluctuation'" << std::endl;
        gnuplotScript << "plot 'rms_energy_fluctuations.dat' with linespoints" << std::endl;
        gnuplotScript << "pause -1" << std::endl;
        gnuplotScript.close();
    } else {
        std::cerr << "Unable to open Gnuplot script file." << std::endl;
    }

    // Print RMS energy fluctuations for each run
    std::cout << "RMS Energy Fluctuations:\n";
    for (int run = 0; run < numRuns; ++run) {
        std::cout << "Run " << (run + 1) << ": " << rmsEnergyFluctuations[run] << std::endl;
    }
}

int main() {
    // Simulation parameters
    int numRuns = 9;                     // Number of runs with different timesteps
    int NRepeats = 10;                     // Number of runs with different timesteps
    std::vector<double> timesteps = {0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};  // Different timesteps for each run

    calculateRMSEnergyFluctuations(NRepeats, numRuns, timesteps);

    return 0;
}
