//? Code for calculating KE per particle for Hard Sphere Model
//? Changed Particle No. & Box Dimension and played arround with the T_step
//? Rewamped the whole code: Used LJ Model instead of Hardsphere, Messed up the froce equation and corrected it, added an initializeSystem Function
//? Corrected the LJ force function ig
//? Added potential energy Function

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <fstream>
#include <ctime>
#include <random>

// Particle structure
struct Particle {
    double x, y, vx, vy, ax, ay;
};

// Function to calculate the inter-particle forces
// Function to calculate the forces and potential energy of the system
double calculateForcesAndEnergy(std::vector<Particle>& particles, double L) {
    double totalEnergy = 0.0;

    // Reset forces and energy
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
    }

    // Calculate forces and potential energy
    for (int i = 0; i < particles.size() - 1; ++i) {
        for (int j = i + 1; j < particles.size(); ++j) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double r = std::sqrt(dx * dx + dy * dy);

            double epsilon = 1.0;  // Depth of the potential well
            double sigma = 1.0;    // Distance at which the potential is zero
            double cutoff = 2.5 * sigma; // Cutoff distance for the Lennard-Jones potential
            double minDistance = 1.0 * sigma; // Minimum distance to prevent particle overlap

            if (r < cutoff && r > minDistance) {
                double r_inv = sigma / r;
                double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                double r_inv12 = r_inv6 * r_inv6;
                double force = 48 * (epsilon/r) * (r_inv12 - 0.5 * r_inv6);

                double fx = force * dx / r;
                double fy = force * dy / r;
                particles[i].ax += fx;
                particles[i].ay += fy;
                particles[j].ax -= fx;
                particles[j].ay -= fy;

                double potential = 4 * epsilon * (r_inv12 - r_inv6);
                totalEnergy += potential;
            }
        }
    }

    return totalEnergy;
}


// Function to update the positions and velocities using the Velocity Verlet method
double updatePositionsAndVelocities(std::vector<Particle>& particles, double dt, double L) {
    double potentialEnergy = calculateForcesAndEnergy(particles, L);

    for (auto& particle : particles) {
        // Update velocities (half step)
        particle.vx += 0.5 * particle.ax * dt;
        particle.vy += 0.5 * particle.ay * dt;

        // Update positions
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;

        // Apply periodic boundary conditions
        particle.x -= std::floor(particle.x / L) * L;
        particle.y -= std::floor(particle.y / L) * L;
    }

    // Update forces and calculate new potential energy
    potentialEnergy = calculateForcesAndEnergy(particles, L);

    for (auto& particle : particles) {
        // Update velocities (half step)
        particle.vx += 0.5 * particle.ax * dt;
        particle.vy += 0.5 * particle.ay * dt;
    }
    return potentialEnergy;
}


// Function to calculate the total kinetic energy per particle of the system
double calculateTotalKineticEnergyPerParticle(const std::vector<Particle>& particles) {
    double totalEnergy = 0.0;
    for (const auto& particle : particles) {
        double vx = particle.vx;
        double vy = particle.vy;
        double kineticEnergy = 0.5 * (vx * vx + vy * vy); // Kinetic energy per particle
        totalEnergy += kineticEnergy;
    }
    return totalEnergy;
}

// Initialize particles
void initializeSystem(std::vector<Particle>& particles, int numParticles, double boxSize, double minSeparation) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (int i = 0; i < numParticles; ++i) {
        bool isOverlap = false;
        do {
            double randomValue1 = dis(gen);
            double randomValue2 = dis(gen);
            particles[i].x = boxSize * randomValue1;
            particles[i].y = boxSize * randomValue2;
            particles[i].vx = 0.0;
            particles[i].vy = 0.0;
            particles[i].ax = 0.0;
            particles[i].ay = 0.0;

            // Check for overlap with previously initialized particles
            for (int j = 0; j < i; ++j) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double separation = std::sqrt(dx * dx + dy * dy);
                if (separation < minSeparation) {
                    isOverlap = true;
                    break;
                }
            }
        } while (isOverlap);
    }
}

int main() {
    // Simulation parameters
    int numParticles = 1000;     // Number of particles         10  20  50 100 500 1000
    double boxSize = 600.0;     // Size of the simulation box   15  30  70 100 350 600 
    double timestep = 0.01;    // Timestep
    int numSteps = 20000;      // Number of simulation steps
    const double minSeparation = 0.7; // Minimum separation between particles (adjust as needed)

    std::vector<Particle> particles(numParticles);
    // Open the output file for writing the data
    initializeSystem(particles, numParticles, boxSize, minSeparation);
    
    std::ofstream outputFile("energy_data.dat");
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open output file." << std::endl;
        return 1;
    }

    // Main simulation loop
    for (int step = 0; step < numSteps; ++step) {

        // Update positions and velocities
        double PotentialEnergyPerParticle=updatePositionsAndVelocities(particles, timestep, boxSize);
        // Calculate and store the total kinetic energy per particle
        double kineticEnergyPerParticle = calculateTotalKineticEnergyPerParticle(particles);
        // Write the data to the output file
        outputFile << step * timestep << " " << kineticEnergyPerParticle << " " << PotentialEnergyPerParticle << std::endl;
    }

    // Close the output file
    outputFile.close();

    // Gnuplot snippet to plot the kinetic energy per particle with time
    std::ofstream gnuplotScript("gnuplot_script.plt");
    if (gnuplotScript.is_open()) {
        gnuplotScript << "set xlabel 'Time'" << std::endl;
        gnuplotScript << "set ylabel 'Energy per Particle'" << std::endl;
        gnuplotScript << "plot 'energy_data.dat' using 1:2 with linespoints title 'Kinetic Energy per Particle', "<< "     'energy_data.dat' using 1:3 with linespoints title 'Potential Energy per Particle'" << std::endl;
        gnuplotScript << "pause -1" << std::endl;
        gnuplotScript.close();
    } else {
        std::cerr << "Unable to open Gnuplot script file." << std::endl;
    }




    // Print a success message
    std::cout << "Simulation completed successfully. Kinetic energy data written to 'energy_data.dat'." << std::endl;

    return 0;
}
