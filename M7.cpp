//added auto params
// Changed saving method to folders based on excecution time
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
#include <algorithm>
#include <string>
#ifdef _WIN32
//compiled on a Windows system
    #include <direct.h>
#else
//compiled on a Unix-based system
    #include <sys/stat.h>
#endif
#include <sys/types.h>
#include <omp.h>
#include "json.hpp" // Include the JSON library, e.g., nlohmann/json

using namespace std;

//params-begin
    //version manger parameters
    int codeVersion=4;
    int currentVersion=0;
    string comment=" Trying for a huge number of particle run, for the displasment and persistence plots, as they are not being exactly logaraithmic, also for the prev run I forgot to change the parameters";
    // Simulation parameters
    double epsilon = 1.0;               // Depth of the potential well
    double sigma = 1.0;                 // Distance at which the potential is zero
    int numParticles = 10;             // !Number of particles          10 100 625 5000 10,000
    double boxSize = 32.0;              // !Size of the simulation box   32 100 250 700  1000
    double temperature = 1.0;           // !Temperature to be fed to MB Distribution
    double cutoff = 1.12246205 * sigma; // Cutoff distance for the Lennard-Jones potential (sigma*2^(1/6)=1.122462048)
    double my_gamma=100.0;              // ! Gamma ;
    //Dynamic Parameters
    double timestep = 0.001;            // !Timestep
    int numSteps = 1000000;              // !Number of simulation steps
    int dataCompression=500;            // If dataCompression-n; The KE & PE at every nth timestep is writen in .dat file
    int period=100;                     // No of time steps to skip before applying next velocity initialization based on MB Distribution
    int NoOfPeriods=1000;                 // No of time velocity intialization is to be applied
    //Vission parameters
    double theta = M_PI/7.0;            // Half of the opening angle of the vision cone (between pi/12 and pi/2)
    double R_0 = 1.5 * sigma;
    double D_R = 0.08;                  // Rotational diffusion coefficient (fixed as such) 0.08
    double D_T= temperature/my_gamma;
    double Pe=(epsilon/temperature)-1;
    double v_0=1.0;                     //(such that Pe= sigma* v_0/D_T=2000
    double Omega = 4.8;                 // Maneuverability strength (60 * D_R)
    //Parellelization parameters
    int numThreads = 8;  // Set the desired number of threads for paralaization

//params-end
//autoset-params-begin
    // double tau=sqrt(sigma*sigma/temperature);
    // double my_gamma=100/tau;
    // double D_R=0.08/tau;
    // double D_T= temperature/my_gamma;
    // double Pe=(epsilon/temperature)-1;
    // double v_0=Pe * D_T / sigma;
    // double Omega=62.5*D_R;
//autoset-params-end


// Particle structure
struct Particle {
    double x, y, vx, vy, ax, ay, phi, PE;
};

// Function to calculate the forces and potential energy of the system
void LJandVissionCone(vector<Particle>& particles, double L, double dt) {
    double totalPotentialEnergy = 0.0; 
    // Reset forces and energy for all particles
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
    }
    // Initialize random number generation for Lambda_i
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> Lambda_i(0.0, 2.0 * D_R);

    // Initialize vectors to store variables related to vision calculations
    std::vector<double> vNc(numParticles, 0.0);
    std::vector<double> vsum_term(numParticles, 0.0);

    // Calculate forces and potential energy between particles (parallelized)
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
        double phi_i = particles[i].phi;
        for (int j = i + 1; j < particles.size(); ++j) {
            double phi_j = particles[j].phi;

            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;

            // Apply minimum image convention to handle particles across periodic boundaries
            if (dx > (L / 2.0)) dx -= L;
            else if (dx < (-L / 2.0)) dx += L;
            if (dy > (L / 2.0)) dy -= L;
            else if (dy < (-L / 2.0)) dy += L;

            double r = sqrt(dx * dx + dy * dy);
            // if (r==0) continue;

            // Lennard-Jones potential and force calculation
            if (r <= cutoff) {
                double r_inv = sigma / r;
                double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                double r_inv12 = r_inv6 * r_inv6;
                double force = -48 * (epsilon / r) * (r_inv12 - 0.5 * r_inv6);

                // Update forces for both particles
                particles[i].ax += force * dx / r;
                particles[i].ay += force * dy / r;
                particles[j].ax -= force * dx / r;
                particles[j].ay -= force * dy / r;

                // Calculate Lennard-Jones potential and accumulate total potential energy
                double potential = 4 * epsilon * (r_inv12 - r_inv6) + epsilon;
                totalPotentialEnergy += potential;
            }

            // Vision Cone calculations
            double phi_ij = atan2(dy, dx);
            double phi_ji = atan2(-dy, -dx);

            if (r <= 4.0 * R_0) {
                // Check if particles are within the vision cone of each other
                if (((dx * cos(phi_i)) + (dy * sin(phi_i))) / r >= cos(theta)) {
                    double delta_phi = phi_ij - phi_i;
                    vsum_term[i] += exp(-r / R_0) * sin(delta_phi);
                    vNc[i] += exp(-r / R_0);
                }
                if (((-dx * cos(phi_j)) + (-dy * sin(phi_j))) / r >= cos(theta)) {
                    double delta_phi = phi_ji - phi_j;
                    vsum_term[j] += exp(-r / R_0) * sin(delta_phi);
                    vNc[j] += exp(-r / R_0);
                }
            }
        }
        

        // Update particle orientation (phi) and angular acceleration (ax, ay)
        double dphi_i;
        if (vNc[i] == 0) continue; // Avoid division by zero
        dphi_i = ((Omega / vNc[i]) * vsum_term[i] + Lambda_i(gen)) * dt; // Calculate angle to update phi_i with 
        phi_i += dphi_i; // Update phi_i
        particles[i].phi = phi_i;
        particles[i].ax += my_gamma * v_0 * cos(particles[i].phi); // Update ax
        particles[i].ay += my_gamma * v_0 * sin(particles[i].phi); // Update ay
    } 
    #pragma omp barrier
    // Calculate and store the average potential energy per particle for a given timestep in the data set of the last particle in the timestep
    particles[numParticles - 1].PE = totalPotentialEnergy / particles.size();
}

// Function to update the positions and velocities using the Velocity Verlet method
void updatePositionsAndVelocities(vector<Particle>& particles, double L , double dt) {

    double k_B=1.0;
    double half_dt = 0.5 * dt;
    double var= sqrt(2.0 * my_gamma * k_B * temperature);
    // Generate random number generator for Gaussian white noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> noiseDist(0.0, 1.0);

    std::vector<std::vector<double>> noise(
        numParticles, std::vector<double>(
            4, 0.0
        )
    );
    
    int i=0;
    for (auto& particle : particles) {
        //Saves Eta and Zeta noise terms in each direction for each particle to be used in thsi and the next for loop
        noise[i][0] = noiseDist(gen);
        noise[i][1] = noiseDist(gen);
        noise[i][2] = noiseDist(gen);
        noise[i][3] = noiseDist(gen);

        // Update velocities (half step)
        // Save the current particle's velocity and acceleration
        double v_xn = particle.vx; double v_yn = particle.vy;
        double a_xn = particle.ax; double a_yn = particle.ay;

        // Update velocities using half step
        particle.vx += half_dt * (a_xn - my_gamma * v_xn) + 0.5 * sqrt(dt) * var * noise[i][0] - 0.125 * dt * dt * my_gamma * (a_xn - my_gamma * v_xn) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][0] + (1.0/sqrt(3.0)) * noise[i][2]);
        particle.vy += half_dt * (a_yn - my_gamma * v_yn) + 0.5 * sqrt(dt) * var * noise[i][1] - 0.125 * dt * dt * my_gamma * (a_yn - my_gamma * v_yn) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][1] + (1.0/sqrt(3.0)) * noise[i][3]);        
        // Update positions
        particle.x += particle.vx * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * noise[i][2];
        particle.y += particle.vy * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * noise[i][3];

        // Apply periodic boundary conditions
        particle.x -= floor(particle.x / L) * L;
        particle.y -= floor(particle.y / L) * L;
        
        i=i+1;

    }
    // Update forces and calculate new potential energy
    LJandVissionCone(particles, L, dt);

    i=0;
    for (auto& particle : particles) {

        double a_x_new = particle.ax; double a_y_new = particle.ay;

        // Update velocities (half step)
        particle.vx += half_dt * (a_x_new - my_gamma * particle.vx) + 0.5 * sqrt(dt) * var * noise[i][0] - 0.125 * dt * dt * my_gamma * (a_x_new - my_gamma * particle.vx) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][0] + (1.0/sqrt(3.0)) * noise[i][2]);
        particle.vy += half_dt * (a_y_new - my_gamma * particle.vy) + 0.5 * sqrt(dt) * var * noise[i][1] - 0.125 * dt * dt * my_gamma * (a_y_new - my_gamma * particle.vy) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][1] + (1.0/sqrt(3.0)) * noise[i][3]);
        i=i+1;
    }
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
        tensor[jsonStep][p][4] = particles[p].phi; // Vission Cone direction 
        if (p==numParticles-1) tensor[jsonStep][p][5] = particles[p].PE; // y-velocity
    }


}



void Notifier(int codeVersion, int currentVersion){
    //Plots the Energy.dat file using plotityy script and feeds it the used parameters which the script uses for naming the created html and png file
    string cd ="python Notifier.py "+ to_string(codeVersion) + " " + to_string(currentVersion);
    system(cd.c_str());
}

void createDirectory(const std::string& dir) {
    // program is being compiled on a Windows system
    #ifdef _WIN32
    //_mkdir function from the direct.h library to create the directory
        _mkdir(dir.c_str());
    // compiled on a Unix-based system
    #else
    // mkdir function from the sys/stat.h library to create the directory. The second argument 0777 sets the permissions of the directory 
        mkdir(dir.c_str(), 0777);
    #endif
}
void save_param(string location){
    // save the above variables as a key value pair in a csv file
    ofstream myfile;
    myfile.open (location+"/params.csv");
    myfile << "comment,"<<comment<<"\n";
    myfile << "epsilon,"<<epsilon<<"\n";
    myfile << "sigma,"<<sigma<<"\n";
    myfile << "numParticles,"<<numParticles<<"\n";
    myfile << "boxSize,"<<boxSize<<"\n";
    myfile << "temperature,"<<temperature<<"\n";
    myfile << "cutoff,"<<cutoff<<"\n";
    myfile << "my_gamma,"<<my_gamma<<"\n";
    myfile << "timestep,"<<timestep<<"\n";
    myfile << "numSteps,"<<numSteps<<"\n";
    myfile << "dataCompression,"<<dataCompression<<"\n";
    myfile << "period,"<<period<<"\n";
    myfile << "NoOfPeriods,"<<NoOfPeriods<<"\n";
    myfile << "theta,"<<theta<<"\n";
    myfile << "R_0,"<<R_0<<"\n";
    myfile << "D_R,"<<D_R<<"\n";
    myfile << "D_T,"<<temperature/my_gamma<<"\n";
    myfile << "Pe,"<<(epsilon/temperature)-1.0<<"\n";
    myfile << "v_0,"<<v_0<<"\n";
    myfile << "Omega,"<<Omega<<"\n";
    myfile.close();
}

string dirlocater(){
    // Get the current time
    time_t now = time(0);

    // Allocate a buffer for the formatted time
    char dt[100];  // Adjust the size based on your formatting needs

    // Format the time using strftime
    strftime(dt, sizeof(dt), "%Y_%m_%d_%H_%M_%S", localtime(&now));

    // Convert ctime output to string and remove newline character
    std::string str_dt(dt);
    str_dt.erase(std::remove(str_dt.end()-1, str_dt.end(), '\n'), str_dt.end());

    // Replace spaces and colons with underscores
    std::replace(str_dt.begin(), str_dt.end(), ' ', '_');
    std::replace(str_dt.begin(), str_dt.end(), ':', '_');
    
    // Create directory
    createDirectory("Runs/"+str_dt);

    return "Runs/"+str_dt;
}

int main() {
    omp_set_num_threads(numThreads);
    // system("python reg.py");
    vector<Particle> particles(numParticles);
    vector<double> velocities;  // Empty velocity vector
    initializeSystem(particles, boxSize);
    vector<vector<vector<double>>> tensor(
        int(floor(double(numSteps)/double(dataCompression))), vector<vector<double>>(
            numParticles, vector<double>(
                6, 0.0
            )
        )
    );
    //Calculating Packing fraction
    double packingFraction = (numParticles*M_PI*(sigma/2)*(sigma/2))/(boxSize*boxSize);
    cout << "Packing Fraction is: " <<packingFraction<<endl;

    string location = dirlocater();
    save_param(location);
    // Main simulation loop
    for (int step = 0; step < numSteps; ++step) {
        // Update positions and velocities
        updatePositionsAndVelocities(particles, boxSize, timestep);

        // Periodicaly initializes the velocity 
        if ( (step % period==0) && (floor(step/period)<NoOfPeriods))
            initializeVelocityAndAcceleration(particles, velocities);

        // Calculate and store the total kinetic energy per particle for every nth data
        if ( step % dataCompression == 0){
            // Write the position data to the output file
            saveParticleData(tensor, particles, step);
        }
        if (fmod((10000.0 * step) / (numSteps * 1.0), 1.0) == 0.0) {
            cout <<"  "<<(100.0 * step) / (numSteps*1.0) << " Percentage completed \r";
        }
    }
    // Open the file in text mode
    std::ofstream file(location+"/particle_positions.json");

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

    // Print a success message
    cout << "Simulation completed successfully. Starting the Plotting and Notifying engines" << endl;


    // Notifier(codeVersion, currentVersion);

    return 0;
}