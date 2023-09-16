//? Added comments
//? 3D-> 2D
//? Added Vission
//? Moved vision to a new function
//? Added gif maker
//? Added break() wherever divission occurs
//? changed round() -> ceil() in particlesPerDimension()
//? Changed break -> continue
//? Added more vaeriables to tensor.json
//? Made global parameters
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
#include "json.hpp" // Include the JSON library, e.g., nlohmann/json

using namespace std;

    //version manger parameters
    int codeVersion=3;
    int currentVersion=0;
    string comment=" Changed Vission Code";
    // Simulation parameters
    double epsilon = 1.0;               // Depth of the potential well
    double sigma = 1.0;                 // Distance at which the potential is zero
    int numParticles = 100;             // !Number of particles          100   625  300
    double boxSize = 50.0;              // !Size of the simulation box   12   250  175
    double temperature = 1.0;           // !Temperature to be fed to MB Distribution
    double my_gamma=100.0;                 // ! Gamma ;
    //Dynamic Parameters
    double timestep = 0.001;            // !Timestep
    int numSteps = 100000;              // !Number of simulation steps
    int dataCompression=500;            // If dataCompression=n; The KE & PE at every nth timestep is writen in .dat file
    int period=100;                     // No of time steps to skip before applying next velocity initialization based on MB Distribution
    int NoOfPeriods=100;                 // No of time velocity intialization is to be applied
    //Vission parameters
    double Omega = 4.8;                 // Maneuverability strength (60 * D_R)
    double D_R = 0.08;                  // Rotational diffusion coefficient (fixed as such) 0.08
    double R_0 = 1.5;                   // Characteristic length      (1.5 * sigma)
    double theta = M_PI/2.0;            // Half of the opening angle of the vision cone (between pi/12 and pi/2)
    double v_0=1.0;                     //(such that Pe= sigma* v_0/D_T=2000

// Particle structure
struct Particle {
    double x, y, vx, vy, ax, ay, zx, zy, ex, ey, phi;
};

// Function to calculate the forces and potential energy of the system
double calculateForcesAndEnergy(vector<Particle>& particles, double L, double dt) {
    double totalPotentialEnergy = 0.0; 
    // Reset forces and energy
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
    }

    double cutoff = 1.12246205 * sigma; // Cutoff distance for the Lennard-Jones potential (sigma*2^(1/6)=1.122462048)

    // Calculate forces and potential energy
    for (int i = 0; i < particles.size() - 1; ++i) {
        for (int j = i + 1; j < particles.size(); ++j) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;

            // Apply minimum image convention
            if (dx > (L / 2.0))dx -= L;
            else if (dx < (-L / 2.0))dx += L;
            if (dy > (L / 2.0))dy -= L;
            else if (dy < (-L / 2.0))dy += L;

            double r = sqrt(dx * dx + dy * dy);
            // if (r==0) continue;

            //Lennard-Jones potential and force calculation 
            if (r <= (cutoff) ) { 
                double r_inv = sigma / r;
                double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                double r_inv12 = r_inv6 * r_inv6;
                double force = -48 * (epsilon / r) * (r_inv12 - 0.5 * r_inv6);
                particles[i].ax += force * dx / r;
                particles[i].ay += force * dy / r;
                particles[j].ax -= force * dx / r;
                particles[j].ay -= force * dy / r;

                double potential = 4 * epsilon * (r_inv12 - r_inv6) + epsilon;
                totalPotentialEnergy += potential;
            }
        }
        particles[i].ax += my_gamma * v_0 * cos(particles[i].phi);
        particles[i].ay += my_gamma * v_0 * sin(particles[i].phi);
    }
    return totalPotentialEnergy/particles.size();
}

void visionCone(vector<Particle>& particles, double L, double dt){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> Lambda_i(0.0,2.0 * D_R);

    // You'll need to incorporate these calculations into your loop over particles
    for (int i = 0; i < particles.size(); ++i) {
        double sum_term = 0.0;
        double Nc_i=0.0;
        double phi_i=particles[i].phi;
        for (int j = 0; j < particles.size(); ++j) {
            if (j != i) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;

                // Apply minimum image convention
                if (dx > (L / 2.0))dx -= L;
                else if (dx < (-L / 2.0))dx += L;
                if (dy > (L / 2.0))dy -= L;
                else if (dy < (-L / 2.0))dy += L;
                double r_ij = sqrt(dx * dx + dy * dy);

                double phi_ij = atan2(dy, dx);
                if ((r_ij <= 4 * R_0) && (((dx * cos(phi_i)) + (dy * sin(phi_i)))/r_ij>= cos(theta))) {
                    double delta_phi = phi_ij - phi_i;
                    sum_term += exp(-r_ij / R_0) * sin(delta_phi);
                    Nc_i += exp(-r_ij / R_0);
                }
            }
        }
        double dphi_i;
        if (Nc_i==0) continue;
        dphi_i = ((Omega / Nc_i) * sum_term + Lambda_i(gen)) * dt; // Update angle phi_i
        phi_i += dphi_i; // Update phi_i
        particles[i].phi=phi_i;
    }
}
// Function to update the positions and velocities using the Velocity Verlet method
double updatePositionsAndVelocities(vector<Particle>& particles, double L , double dt) {

    double k_B=1.0;
    double half_dt = 0.5 * dt;
    double var= sqrt(2.0 * my_gamma * k_B * temperature);
    // Generate random number generator for Gaussian white noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> noiseDist(0.0, 1.0);
    

    for (auto& particle : particles) {
        //Saves Eta and Zeta noise terms in each direction for each particle to be used in thsi and the next for loop
        particle.ex = noiseDist(gen);
        particle.ey = noiseDist(gen);
        particle.zx = noiseDist(gen);
        particle.zy = noiseDist(gen);

        // Update velocities (half step)
        // Save the current particle's velocity and acceleration
        double v_xn = particle.vx; double v_yn = particle.vy;
        double a_xn = particle.ax; double a_yn = particle.ay;

        // Update velocities using half step
        particle.vx += half_dt * (a_xn - my_gamma * v_xn) + 0.5 * sqrt(dt) * var * particle.ex - 0.125 * dt * dt * my_gamma * (a_xn - my_gamma * v_xn) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * particle.ex + (1.0/sqrt(3.0)) * particle.zx);
        particle.vy += half_dt * (a_yn - my_gamma * v_yn) + 0.5 * sqrt(dt) * var * particle.ey - 0.125 * dt * dt * my_gamma * (a_yn - my_gamma * v_yn) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * particle.ey + (1.0/sqrt(3.0)) * particle.zy);        
        // Update positions
        particle.x += particle.vx * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * particle.zx;
        particle.y += particle.vy * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * particle.zy;

        // Apply periodic boundary conditions
        particle.x -= floor(particle.x / L) * L;
        particle.y -= floor(particle.y / L) * L;
        

    }
    visionCone(particles, L, dt);
    // Update forces and calculate new potential energy
    double potentialEnergyPerParticle = calculateForcesAndEnergy(particles, L, dt);

    for (auto& particle : particles) {

        double a_x_new = particle.ax; double a_y_new = particle.ay;

        // Update velocities (half step)
        particle.vx += half_dt * (a_x_new - my_gamma * particle.vx) + 0.5 * sqrt(dt) * var * particle.ex - 0.125 * dt * dt * my_gamma * (a_x_new - my_gamma * particle.vx) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * particle.ex + (1.0/sqrt(3.0)) * particle.zx);
        particle.vy += half_dt * (a_y_new - my_gamma * particle.vy) + 0.5 * sqrt(dt) * var * particle.ey - 0.125 * dt * dt * my_gamma * (a_y_new - my_gamma * particle.vy) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * particle.ey + (1.0/sqrt(3.0)) * particle.zy);
    }
    return potentialEnergyPerParticle;
}

// Function to calculate the total kinetic energy per particle of the system
double calculateTotalKineticEnergyPerParticle(const vector<Particle>& particles) {
    double totalKineticEnergy = 0.0;
    for (const auto& particle : particles) {
        double vx = particle.vx;
        double vy = particle.vy;
        double kineticEnergy = 0.5 * (vx * vx + vy * vy); // Kinetic energy per particle
        totalKineticEnergy += kineticEnergy;
    }
    return totalKineticEnergy/particles.size();
}

//Snippet to initialize particle possition in a cubic lattice
void initializeSystem(vector<Particle>& particles, double L) {
    random_device rd;
    mt19937 gen(rd());
    //we take -Pi to Pi instead of 0 to 2Pi cuz atan() gives result in the former way
    std::uniform_real_distribution<double> uniform(-M_PI,M_PI);
    // Calculate the number of particles per dimension in the cubic grid
    int particlesPerDimension = ceil(sqrt(numParticles));

    // Calculate the spacing between particles based on the separation distance
    double spacing = L/double(particlesPerDimension);

    // Resize the vector to hold the desired number of particles
    particles.resize(numParticles);

    int index = 0;

    // Loop over each dimension of the cubic grid (x, y, z)
    for (int x = 0; x < particlesPerDimension; ++x) {
        for (int y = 0; y < particlesPerDimension; ++y) {
            // Check if all particles have been placed, and break the loop if so
            if (index >= numParticles) break;

            // Calculate the coordinates of the particle within the cubic grid, The "+ 0.5" is used to center the particle in each grid cell
            particles[index].x = (x + 0.5) * spacing;
            particles[index].y = (y + 0.5) * spacing;
            //Randomnly assigns each particle a vission cone direction
            particles[index].phi = uniform(gen);

            index++;
        }
    }
}

void initializeVelocityAndAcceleration(vector<Particle>& particles, vector<double>& velocities) {
    int numParticles = particles.size();

    // Set random seed
    random_device rd;
    mt19937 gen(rd());

    // Boltzmann factor times temperature
    const double k_T = temperature;

    // Setup the Maxwell distribution, i.e., gamma distribution with alpha = 1.0
    gamma_distribution<double> maxwell(1.0, k_T);

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
        double v = sqrt(particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
        if (v==0) continue ;
        particles[i].vx = speed * (particles[i].vx / v);
        particles[i].vy = speed * (particles[i].vy / v);
    }
    //REVIEW - Why is this here
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
    }
}

void saveParticleData( std::vector<std::vector<std::vector<double>>>& tensor, std::vector<Particle> particles, int step) {
    // Get the number of time steps, particles, and dimensions
    int numParticles = static_cast<int>(particles.size());

    int jsonStep=int(floor(double(step)/double(dataCompression)));
    // Fill the tensor with particle positions
    for (int p = 0; p < numParticles; ++p) {
        tensor[jsonStep][p][0] = particles[p].x; // x-coordinate
        tensor[jsonStep][p][1] = particles[p].y; // y-coordinate
        tensor[jsonStep][p][2] = particles[p].vx; // x-velocity
        tensor[jsonStep][p][3] = particles[p].vy; // y-velocity
        tensor[jsonStep][p][4] = particles[p].phi; // y-velocity
    }

    // Open the file in text mode
        std::ofstream file("particle_positions.json");

        if (file.is_open()) {
            // Serialize the tensor to JSON
            nlohmann::json jsonData(tensor);

            // Write the JSON data to the file
            file << jsonData.dump(4); // The argument sets the indentation for readability

            // Close the file
            file.close();
        } else {
            std::cerr << "Unable to open the JSON file for writing." << std::endl;
        }
}

int versionUpdater(){
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
    //Plots the Energy.dat file using plotityy script and feeds it the used parameters which the script uses for naming the created html and png file
    string cd ="python plotityy.py "+ to_string(codeVersion) + " " + to_string(currentVersion)+" " + to_string(numParticles) +" "+to_string(int(boxSize))+" "+to_string(int(temperature))+" "+to_string(timestep)+" "+ to_string(numSteps);
    system(cd.c_str());
    //Creates a file name based on used parameters
    string filename="M6-V"+ to_string(codeVersion) + "." + to_string(currentVersion)+"_nP-"+ to_string(numParticles) +"_Bs-"+to_string(int(boxSize))+"_T-"+to_string(int(temperature))+"_ls-"+to_string(timestep)+"_ns-"+ to_string(numSteps);
    //Finds the file based on the file name and auto-sends it to my mail
    string ced= "python sendmail.py \"Graph Generated\" -b \"" + comment +"\" -a "+filename + ".png,"+ filename+".html" ;
    system(ced.c_str());
}

int main() {
    vector<Particle> particles(numParticles);
    vector<double> velocities;  // Empty velocity vector
    initializeSystem(particles, boxSize);
    std::vector<std::vector<std::vector<double>>> tensor(
        int(floor(double(numSteps)/double(dataCompression))), std::vector<std::vector<double>>(
            numParticles, std::vector<double>(
                5, 0.0
            )
        )
    );
    //Calculating Packing fraction
    double packingFraction = (numParticles*M_PI*(sigma/2)*(sigma/2))/(boxSize*boxSize);
    cout << "Packing Fraction is: " <<packingFraction<<endl;
    
    // Open the output file for writing the data
    ofstream outputFile("energy_data.dat");
    if (!outputFile.is_open()) {
        cerr << "Unable to open output file." << endl;
        return 1;
    }
    // Main simulation loop
    for (int step = 0; step < numSteps; ++step) {
        // Update positions and velocities
        double PotentialEnergyPerParticle = updatePositionsAndVelocities(particles, boxSize, timestep);
        

        // Periodicaly initializes the velocity 
        if ( (step % period==0) && (floor(step/period)<NoOfPeriods))
            initializeVelocityAndAcceleration(particles, velocities);

        // Calculate and store the total kinetic energy per particle for every nth data
        if ( step % dataCompression == 0){
            // Write the position data to the output file
            saveParticleData(tensor, particles, step);
            
            double kineticEnergyPerParticle = calculateTotalKineticEnergyPerParticle(particles);
            // Write the Energy data to the output file
            outputFile << step * timestep << " " << kineticEnergyPerParticle << " " << PotentialEnergyPerParticle << " " << (kineticEnergyPerParticle+PotentialEnergyPerParticle) << endl;
        }
        if (fmod((1000.0 * step) / (numSteps * 1.0), 1.0) == 0.0) {
            cout <<"  "<<(100.0 * step) / (numSteps*1.0) << " Percentage completed \r";
        }
    }
    // Close the output file
    outputFile.close();

    // Print a success message
    cout << "Simulation completed successfully." << endl;


    // finalPositionRecorder(particles);
    currentVersion=versionUpdater();
    plotterCumMailer(codeVersion, currentVersion, numParticles, boxSize, temperature, timestep, numSteps, comment);

    return 0;
}