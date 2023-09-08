//? Changed the code for MB distribution from velocity to energy (Fixed the KE getting fixed at random energy issue)
//? Added automatic plot version updater
//? Moved the fInal possition recorder to an external function
//? Added white Noice
//? Added versionUpdater and plotterCumMailer as function
//? Added Euler–Maruyama scheme
//? Fixed r_inv6 value
//? Changed the min sep value to boxSize/double(particlesPerDimension), therby fixed the blowups at high PF
//? Properly assigned noise
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <fstream>
#include <ctime>
#include <random>
#include <string>

using namespace std;

// Particle structure
struct Particle {
    double x, y, z, vx, vy, vz, ax, ay, az, zx, zy, zz, ex, ey, ez;
};

// Function to calculate the forces and potential energy of the system
double calculateForcesAndEnergy(vector<Particle>& particles, double L, double sigma, double epsilon, double temperature, double gamma, double dt) {
    double totalPotentialEnergy = 0.0;

    // Reset forces and energy
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
        particle.az = 0.0;
    }

    double cutoff = 2.5 * sigma; // Cutoff distance for the Lennard-Jones potential (sigma*2^(1/6)=1.122462048)

    
    // Calculate forces and potential energy
    for (int i = 0; i < particles.size() - 1; ++i) {
        for (int j = i + 1; j < particles.size(); ++j) {
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            // Apply minimum image convention
            if (dx > L / 2)dx -= L;
            else if (dx < -L / 2)dx += L;

            if (dy > L / 2)dy -= L;
            else if (dy < -L / 2)dy += L;

            if (dz > L / 2)dz -= L;
            else if (dz < -L / 2)dz += L;

            double r = sqrt(dx * dx + dy * dy + dz * dz);

            //Lennard-Jones potential and force calculation 
            if (r <= cutoff ) { 
                double r_inv = sigma / r;
                double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                double r_inv12 = r_inv6 * r_inv6;
                double force = -48 * (epsilon / r) * (r_inv12 - 0.5 * r_inv6);
                // Add Gaussian white noise to forces

                particles[i].ax -= force * dx / r;
                particles[i].ay -= force * dy / r;
                particles[i].az -= force * dz / r;
                particles[j].ax += force * dx / r;
                particles[j].ay += force * dy / r;
                particles[j].az += force * dz / r;


                double potential = 4 * epsilon * (r_inv12 - r_inv6);
                totalPotentialEnergy += potential;
            }
        }
    }

    return totalPotentialEnergy/particles.size();
}

// Function to update the positions and velocities using the Velocity Verlet method
double updatePositionsAndVelocities(vector<Particle>& particles, double dt, double L,double sigma, double epsilon, double temperature, double gamma) {

    double k_B=1.0;
    double half_dt = 0.5 * dt;
    double var= sqrt(2.0 * gamma * k_B * temperature);
    // Generate random number generator for Gaussian white noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> noiseDist(0.0, 1.0);
    

    for (auto& particle : particles) {
        particle.ex = noiseDist(gen);
        particle.ey = noiseDist(gen);
        particle.ez = noiseDist(gen);
        particle.zx = noiseDist(gen);
        particle.zy = noiseDist(gen);
        particle.zz = noiseDist(gen);

        // Update velocities (half step)
        // Save the current particle's velocity and acceleration
        double v_xn = particle.vx; double v_yn = particle.vy; double v_zn = particle.vz;
        double a_xn = particle.ax; double a_yn = particle.ay; double a_zn = particle.az;

        // Update velocities using half step
        particle.vx += half_dt * (a_xn - gamma * v_xn) + 0.5 * sqrt(dt) * var * particle.ex - 0.125 * dt * dt * gamma * (a_xn - gamma * v_xn) - 0.25 * sqrt(dt) * dt * gamma * var * (0.5 * particle.ex + (1.0/sqrt(3.0)) * particle.zx);
        particle.vy += half_dt * (a_yn - gamma * v_yn) + 0.5 * sqrt(dt) * var * particle.ey - 0.125 * dt * dt * gamma * (a_yn - gamma * v_yn) - 0.25 * sqrt(dt) * dt * gamma * var * (0.5 * particle.ey + (1.0/sqrt(3.0)) * particle.zy);
        particle.vz += half_dt * (a_zn - gamma * v_zn) + 0.5 * sqrt(dt) * var * particle.ez - 0.125 * dt * dt * gamma * (a_zn - gamma * v_zn) - 0.25 * sqrt(dt) * dt * gamma * var * (0.5 * particle.ez + (1.0/sqrt(3.0)) * particle.zz);
        
        // Update positions
        particle.x += particle.vx * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * particle.zx;
        particle.y += particle.vy * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * particle.zy;
        particle.z += particle.vz * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * particle.zz;

        // Apply periodic boundary conditions
        particle.x -= floor(particle.x / L) * L;
        particle.y -= floor(particle.y / L) * L;
        particle.z -= floor(particle.z / L) * L;
        

    }

    // Update forces and calculate new potential energy
    double potentialEnergyPerParticle = calculateForcesAndEnergy(particles, L, sigma, epsilon, temperature, gamma, dt);

    for (auto& particle : particles) {

        double a_x_new = particle.ax;double a_y_new = particle.ay;double a_z_new = particle.az;

        // Update velocities (half step)
        particle.vx += half_dt * (a_x_new - gamma * particle.vx) + 0.5 * sqrt(dt) * var * particle.ex - 0.125 * dt * dt * gamma * (a_x_new - gamma * particle.vx) - 0.25 * sqrt(dt) * dt * gamma * var * (0.5 * particle.ex + (1.0/sqrt(3.0)) * particle.zx);
        particle.vy += half_dt * (a_y_new - gamma * particle.vy) + 0.5 * sqrt(dt) * var * particle.ey - 0.125 * dt * dt * gamma * (a_y_new - gamma * particle.vy) - 0.25 * sqrt(dt) * dt * gamma * var * (0.5 * particle.ey + (1.0/sqrt(3.0)) * particle.zy);
        particle.vz += half_dt * (a_z_new - gamma * particle.vz) + 0.5 * sqrt(dt) * var * particle.ez - 0.125 * dt * dt * gamma * (a_z_new - gamma * particle.vz) - 0.25 * sqrt(dt) * dt * gamma * var * (0.5 * particle.ez + (1.0/sqrt(3.0)) * particle.zz);
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

//Snippet to initialize particle possition in a cubic lattice
void initializeSystem(vector<Particle>& particles, int numParticles, double boxSize) {
    // Calculate the number of particles per dimension in the cubic grid
    int particlesPerDimension = round(cbrt(numParticles));

    // Calculate the spacing between particles based on the separation distance
    double spacing = boxSize/double(particlesPerDimension);

    // Resize the vector to hold the desired number of particles
    particles.resize(numParticles);

    int index = 0;

    // Loop over each dimension of the cubic grid (x, y, z)
    for (int x = 0; x < particlesPerDimension; ++x) {
        for (int y = 0; y < particlesPerDimension; ++y) {
            for (int z = 0; z < particlesPerDimension; ++z) {
                // Check if all particles have been placed, and break the loop if so
                if (index >= numParticles) break;

                // Calculate the coordinates of the particle within the cubic grid, The "+ 0.5" is used to center the particle in each grid cell
                particles[index].x = (x + 0.5) * spacing;
                particles[index].y = (y + 0.5) * spacing;
                particles[index].z = (z + 0.5) * spacing;

                index++;
            }
        }
    }
}

void initializeVelocityAndAcceleration(vector<Particle>& particles, vector<double>& velocities, double temperature) {
    int numParticles = particles.size();

    // Set random seed
    random_device rd;
    mt19937 gen(rd());

    // Boltzmann factor times temperature
    const double k_T = temperature;

    // Setup the Maxwell distribution, i.e., gamma distribution with alpha = 3/2
    gamma_distribution<double> maxwell(3.0 / 2.0, k_T);

    // If velocities vector is empty, generate new velocities
    if (velocities.empty()) {
        velocities.resize(numParticles);
        for (int i = 0; i < numParticles; ++i) {
            // Generate velocity from the Maxwell-Boltzmann distribution
            velocities[i] = sqrt(2 * maxwell(gen));
        }
    }

    // Assign velocities to particles
    for (int i = 0; i < numParticles; ++i) {
        double speed = velocities[i];
        double v = sqrt(particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy + particles[i].vz * particles[i].vz);

        particles[i].vx = speed * (particles[i].vx / v);
        particles[i].vy = speed * (particles[i].vy / v);
        particles[i].vz = speed * (particles[i].vz / v);
    }
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
        particle.az = 0.0;
    }
}

void finalPositionRecorder(vector<Particle> particles){

        // Record final positions of particles
    ofstream finalPositionsFile("final_positions.dat");
    if (finalPositionsFile.is_open()) {
        for (const auto& particle : particles) {
            finalPositionsFile << particle.x << " " << particle.y << " " << particle.z << endl;
        }
        finalPositionsFile.close();
        cout << "Final positions of particles written to 'final_positions.dat'." << endl;
    } else {
        cerr << "Unable to open final positions file." << endl;
    }
}

int versionUpdater(int currentVersion){
    int currentPlotVersion = currentVersion;
    std::ifstream versionNumber("version.dat"); // Open the file for reading
    if (versionNumber.is_open()) {
        versionNumber >> currentPlotVersion; // Read the current integer from the file
        versionNumber.close();

        // Increment the integer by one
        int newVersion = currentPlotVersion + 1;

        std::ofstream versionNumberOut("version.dat"); // Open the file for writing
        if (versionNumberOut.is_open()) {
            versionNumberOut << newVersion; // Write the new integer to the file
            versionNumberOut.close();
            std::cout << "Updated version: " << newVersion << std::endl;
        } else {
            std::cerr << "Unable to open version file for writing." << std::endl;
        }
    } else {
        std::cerr << "Unable to open version file for reading." << std::endl;
    }
    return currentPlotVersion;
}

void plotterCumMailer(int codeVersion, int currentVersion, int (numParticles), double boxSize, double temperature, double timestep, int numSteps, string comment){

    string cd ="python plotityy.py "+ to_string(codeVersion) + " " + to_string(currentVersion)+" " + to_string(numParticles) +" "+to_string(int(boxSize))+" "+to_string(int(temperature))+" "+to_string(timestep)+" "+ to_string(numSteps);
    system(cd.c_str());
    cout << cd.c_str();
    string filename="M5-V"+ to_string(codeVersion) + "." + to_string(currentVersion)+"_nP-"+ to_string(numParticles) +"_Bs-"+to_string(int(boxSize))+"_T-"+to_string(int(temperature))+"_ls-"+to_string(timestep)+"_ns-"+ to_string(numSteps);
    string ced= "python sendmail.py \"Graph Generated\" -b \"" + comment +"\" -a "+filename + ".png,"+ filename+".html" ;
    system(ced.c_str());
}

int main() {
    //version manger parameters
    int codeVersion=5;
    int currentVersion=0;
    string comment=" ";
    // Simulation parameters
    double epsilon = 1.0;               // Depth of the potential well
    double sigma = 1.0;                 // Distance at which the potential is zero
    int numParticles = 100;             // !Number of particles          10 100 1000
    double boxSize = 5.0;               // !Size of the simulation box   2  5   10
    double temperature = 1.0;           // !Temperature to be fed to MB Distribution
    double gamma=1.0;                  // ! Gamma 
    double timestep = 0.0001;           // !Timestep
    int numSteps = 3000000;              // !Number of simulation steps
    int dataCompression=100;            // If dataCompression=n; The KE & PE at every nth timestep is writen in .dat file
    int period=100;                     // No of time steps to skip before applying next velocity initialization based on MB Distribution
    int NoOfPeriods=3000;               // No of time velocity intialization is to be applied
     // double minSeparation = 1.0;         // Minimum separation between particles (adjust as needed)

    double packingFraction = (numParticles*(4.0/3.0)*M_PI*(sigma/2)*(sigma/2)*(sigma/2))/(boxSize*boxSize*boxSize);
    cout << "Packing Fraction is: " <<packingFraction<<endl;

    vector<Particle> particles(numParticles);
    vector<double> velocities;  // Empty velocity vector
    initializeSystem(particles, numParticles, boxSize);
    
    // Open the output file for writing the data
    ofstream outputFile("energy_data.dat");
    if (!outputFile.is_open()) {
        cerr << "Unable to open output file." << endl;
        return 1;
    }
    // Main simulation loop
    for (int step = 0; step < numSteps; ++step) {
        // Update positions and velocities
        double PotentialEnergyPerParticle = updatePositionsAndVelocities(particles, timestep, boxSize, sigma, epsilon, temperature, gamma);

        // Periodicaly initializes the velocity
        if ( step % period==0 && floor(step/period)<NoOfPeriods)
            initializeVelocityAndAcceleration(particles, velocities, temperature);

        // Calculate and store the total kinetic energy per particle for every nth data
        if ( step % dataCompression == 0){
            double kineticEnergyPerParticle = calculateTotalKineticEnergyPerParticle(particles);
            // Write the data to the output file
            outputFile << step * timestep << " " << kineticEnergyPerParticle << " " << PotentialEnergyPerParticle << " " << (kineticEnergyPerParticle+PotentialEnergyPerParticle) << endl;
        }
        if (fmod((1000.0 * step) / (numSteps * 1.0), 1.0) == 0.0) {
            cout <<"  "<<(100.0 * step) / (numSteps*1.0) << " Percentage completed \r";
        }
    }
    // Close the output file
    outputFile.close();
    // Print a success message
    cout << "Simulation completed successfully. Kinetic energy data written to 'energy_data.dat'." << endl;


    finalPositionRecorder(particles);
    currentVersion=versionUpdater(currentVersion);
    plotterCumMailer(codeVersion, currentVersion, numParticles, boxSize, temperature, timestep, numSteps, comment);

    return 0;
}