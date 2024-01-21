//? Added minimum image convention
//? Initialize system updated from random initialization to cubic initialization
//? Also ploted wrong total energy plot
//? Added data compression to increase speed while plotting
//? Updated from 2D to 3D
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <fstream>
#include <ctime>
#include <random>

using namespace std;

// Particle structure
struct Particle {
    double x, y, z, vx, vy, vz, ax, ay, az;
};

// Function to calculate the inter-particle forces
// Function to calculate the forces and potential energy of the system
double calculateForcesAndEnergy(vector<Particle>& particles, double L) {
    double totalPotentialEnergy = 0.0;

    // Reset forces and energy
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
        particle.az = 0.0;
    }

    
    double epsilon = 1.0;  // Depth of the potential well
    double sigma = 1.0;    // Distance at which the potential is zero
    double cutoff = 2.5 * sigma; // Cutoff distance for the Lennard-Jones potential
    double minDistance = 1.0 * sigma; // Minimum distance to prevent particle overlap


    // Calculate forces and potential energy
    for (int i = 0; i < particles.size() - 1; ++i) {
        for (int j = i + 1; j < particles.size(); ++j) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;

            // Apply minimum image convention
            if (dx > L / 2)
                dx -= L;
            else if (dx < -L / 2)
                dx += L;

            if (dy > L / 2)
                dy -= L;
            else if (dy < -L / 2)
                dy += L;

            if (dz > L / 2)
                dz -= L;
            else if (dz < -L / 2)
                dz += L;

            double r = sqrt(dx * dx + dy * dy + dz * dz);

            if (r < cutoff && r > minDistance) {
                double r_inv = sigma / r;
                double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                double r_inv12 = r_inv6 * r_inv6;
                double force = 48 * (epsilon / r) * (r_inv12 - 0.5 * r_inv6);

                double fx = force * dx / r;
                double fy = force * dy / r;
                double fz = force * dz / r;
                particles[i].ax += fx;
                particles[i].ay += fy;
                particles[i].az += fz;
                particles[j].ax -= fx;
                particles[j].ay -= fy;
                particles[j].az -= fz;

                double potential = 4 * epsilon * (r_inv12 - r_inv6);
                totalPotentialEnergy += potential;
            }
        }
    }

    return totalPotentialEnergy/particles.size();
}

// Function to update the positions and velocities using the Velocity Verlet method
double updatePositionsAndVelocities(vector<Particle>& particles, double dt, double L) {
    double potentialEnergyPerParticle = calculateForcesAndEnergy(particles, L);

    for (auto& particle : particles) {
        // Update velocities (half step)
        particle.vx += 0.5 * particle.ax * dt;
        particle.vy += 0.5 * particle.ay * dt;
        particle.vz += 0.5 * particle.az * dt;

        // Update positions
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
        particle.z += particle.vz * dt;

        // Apply periodic boundary conditions
        particle.x -= floor(particle.x / L) * L;
        particle.y -= floor(particle.y / L) * L;
        particle.z -= floor(particle.z / L) * L;
    }

    // Update forces and calculate new potential energy
    potentialEnergyPerParticle = calculateForcesAndEnergy(particles, L);

    for (auto& particle : particles) {
        // Update velocities (half step)
        particle.vx += 0.5 * particle.ax * dt;
        particle.vy += 0.5 * particle.ay * dt;
        particle.vz += 0.5 * particle.az * dt;
    }
    return potentialEnergyPerParticle;
}

// Function to calculate the total kinetic energy per particle of the system
double calculateTotalKineticEnergyPerParticle(const vector<Particle>& particles) {
    double totalKineticEnergy = 0.0;
    for (const auto& particle : particles) {
        double vx = particle.vx;
        double vy = particle.vy;
        double vz = particle.vz;
        double kineticEnergy = 0.5 * (vx * vx + vy * vy + vz * vz); // Kinetic energy per particle
        totalKineticEnergy += kineticEnergy;
    }
    return totalKineticEnergy/particles.size();
}

void initializeSystem(vector<Particle>& particles, int numParticles, double boxSize, double separation) {
    int particlesPerDimension = round(cbrt(numParticles));
    double spacing = separation;

    particles.resize(numParticles);

    int index = 0;
    for (int x = 0; x < particlesPerDimension; ++x) {
        for (int y = 0; y < particlesPerDimension; ++y) {
            for (int z = 0; z < particlesPerDimension; ++z) {
                if (index >= numParticles)
                    break;

                particles[index].x = (x + 0.5) * spacing;
                particles[index].y = (y + 0.5) * spacing;
                particles[index].z = (z + 0.5) * spacing;

                index++;
            }
        }
    }
}

int main() {
    // Simulation parameters
    int numParticles = 100;     // Number of particles          10 100
    double boxSize = 5.0;     // Size of the simulation box    2  5
    double timestep = 0.0001;    // Timestep
    int numSteps = 1000000;      // Number of simulation steps
    int dataCompression=1;
    const double minSeparation = 1.0; // Minimum separation between particles (adjust as needed)

    vector<Particle> particles(numParticles);
    // Open the output file for writing the data
    initializeSystem(particles, numParticles, boxSize, minSeparation);

    ofstream outputFile("energy_data.dat");
    if (!outputFile.is_open()) {
        cerr << "Unable to open output file." << endl;
        return 1;
    }

    // Main simulation loop
    for (int step = 0; step < numSteps; ++step) {
        // Update positions and velocities
        double PotentialEnergyPerParticle = updatePositionsAndVelocities(particles, timestep, boxSize);
        // Calculate and store the total kinetic energy per particle
        // Write the data to the output file
        if ( step % dataCompression == 0){
            double kineticEnergyPerParticle = calculateTotalKineticEnergyPerParticle(particles);
            outputFile << step * timestep << " " << kineticEnergyPerParticle << " " << PotentialEnergyPerParticle << " " << (kineticEnergyPerParticle-PotentialEnergyPerParticle) << endl;
        }
    }

    // Close the output file
    outputFile.close();

    // Gnuplot snippet to plot the kinetic energy per particle with time
    ofstream gnuplotScript("gnuplot_script.plt");
    if (gnuplotScript.is_open()) {
        gnuplotScript << "set xlabel 'Time'" << endl;
        gnuplotScript << "set ylabel 'Energy per Particle'" << endl;
        gnuplotScript << "plot 'energy_data.dat' using 1:2 with linespoints title 'Kinetic Energy per Particle', "
        << "'energy_data.dat' using 1:3 with linespoints title 'Potential Energy per Particle',"
        << "'energy_data.dat' using 1:4 with linespoints title 'Total Energy per Particle'" << endl;
        gnuplotScript << "pause -1" << endl;
        gnuplotScript.close();
    } else {
        cerr << "Unable to open Gnuplot script file." << endl;
    }

    // Print a success message
    cout << "Simulation completed successfully. Kinetic energy data written to 'energy_data.dat'." << endl;

    return 0;
}
