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
    int codeVersion = 0;
    int currentVersion = 0;
    string comment="add comment";
    // Simulation parameters
    double particleRatio= 1.0;          //Percentage of particle1 in total
    double epsilon = 1.0;               // Depth of the potential well
    double sigma = 1.0;                 // Distance at which the potential is zero
    int numParticles = 625;             // *Number of particles
    // double boxSize = 100.0;          // Size of the simulation box
    double packingFraction = 0.0785;
    double boxSize = sqrt((numParticles*M_PI*(sigma/2)*(sigma/2))/packingFraction);
    double cutoff = 1.12246205 * sigma; // Cutoff distance for the Lennard-Jones potential (sigma*2^(1/6)=1.122462048)
    //Dynamic Parameters
    // double timestep = 0.001;         // Timestep
    int numSteps = 1000000;             // *Number of simulation steps
    int period = 100;                   // No of time steps to skip before applying next velocity initialization based on MB Distribution
    int NoOfPeriods=numSteps/1000;      // No of time velocity intialization is to be applied
    int dataCompression = numSteps/1000;// If dataCompression-n; The KE & PE at every nth timestep is writen in .dat file
    //Vission parameters
    double theta = 36.0*(M_PI/180.0);   // *Half of the opening angle of the vision cone (between pi/12 and pi/2)
    double R_0 = 1.5 * sigma;
    // double Pe=(epsilon/temperature)-1;
//params-end
//autoset-params-begin
    double Pe=200.0;
    double temperature=epsilon/(Pe+1.0);
    double tau=sqrt(sigma*sigma/temperature);
    double timestep=0.001*tau;
    double my_gamma=100.0/tau;
    double D_R=0.08/tau;
    double D_T= temperature/my_gamma;
    double v_0=Pe * D_T / sigma;
    double Omega=62.5*D_R;
//autoset-params-end


// Particle structure
struct Particle {
    int particleID;
    double m ,x, y, vx, vy, ax, ay, phi, PE;
};

// Function to calculate the forces and potential energy of the system
void LJandVissionCone(vector<Particle>& particles, vector<double>& vNc, vector<double>& vsum_term, double L) {
    double totalPotentialEnergy = 0.0; 
    // Reset forces and energy for all particles
    for (auto& particle : particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;
    }
    // set vNc and vsum_term to zero for all particles
    vNc.assign(particles.size(), 0.0);
    vsum_term.assign(particles.size(), 0.0);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> noiseDist(0.0, 1.0);
    normal_distribution<double> Lambda_i(0.0, sqrt(2.0 * D_R));

    // Create cell list
    double cellSize = 4.0 * R_0; // cell size should be at least the cutoff distance
    int numCellsPerSide = ceil(L / cellSize);
    vector<vector<int>> cells(numCellsPerSide * numCellsPerSide);

    // Assign particles to cells
    for (int i = 0; i < particles.size(); ++i) {
        int cellX = particles[i].x / cellSize;
        int cellY = particles[i].y / cellSize;
        int cellIndex = cellX + cellY * numCellsPerSide;
        cells[cellIndex].push_back(i);
    }

    // Calculate forces and potential energy between particles (parallelized)
    for (int cellIndex = 0; cellIndex < cells.size(); ++cellIndex) {
        for (int iIndex = 0; iIndex < cells[cellIndex].size(); ++iIndex) {
            int i = cells[cellIndex][iIndex];
            double phi_i = particles[i].phi;
            double dphi_i = 0.0;
            // Check particles in the same cell and neighboring cells
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    int neighborCellX = (cellIndex % numCellsPerSide + dx + numCellsPerSide) % numCellsPerSide;
                    int neighborCellY = (cellIndex / numCellsPerSide + dy + numCellsPerSide) % numCellsPerSide;
                    int neighborCellIndex = neighborCellX + neighborCellY * numCellsPerSide;

                    for (int jIndex = 0; jIndex < cells[neighborCellIndex].size(); ++jIndex) {
                        int j = cells[neighborCellIndex][jIndex];
                        if (i == j) continue; // Skip self-interaction

                        double dx = particles[j].x - particles[i].x;
                        double dy = particles[j].y - particles[i].y;

                        // Apply minimum image convention to handle particles across periodic boundaries
                        dx -= round(dx / L) * L;
                        dy -= round(dy / L) * L;
                        double r = sqrt(dx * dx + dy * dy);

                        // Lennard-Jones potential and force calculation
                        if (r <= cutoff) {
                            double r_inv = sigma / r;
                            double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                            double r_inv12 = r_inv6 * r_inv6;
                            double force = 48 * (epsilon / r) * (r_inv12 - 0.5 * r_inv6);

                            // Update forces for both particles
                            particles[i].ax -= force * dx / r;
                            particles[i].ay -= force * dy / r;

                            // Calculate Lennard-Jones potential and accumulate total potential energy
                            double potential = 4 * epsilon * (r_inv12 - r_inv6) + epsilon;
                            totalPotentialEnergy += potential;
                        }

                        // Vision Cone calculations
                        if (r <= 4.0 * R_0 && ((dx * cos(phi_i)) + (dy * sin(phi_i))) / r >= cos(theta)) {
                            double phi_ij = atan2(dy, dx);
                            double delta_phi = phi_ij - phi_i;
                            vsum_term[i] += exp(-r / R_0) * sin(delta_phi);
                            vNc[i] += exp(-r / R_0);
                        }
                    }
                }
            }
            // Update particle orientation (phi) and angular acceleration (ax, ay)
            if (vNc[i] == 0.0) {dphi_i=Lambda_i(gen)* sqrt(timestep);} // Avoid division by zero
            else {dphi_i = (Omega / vNc[i]) * vsum_term[i] * timestep + Lambda_i(gen)* sqrt(timestep);} // Calculate angle to update phi_i with i
            particles[i].phi += dphi_i;
            // if id = 1, then update the acceleration
            if (particles[i].particleID == 1) {
                particles[i].ax += my_gamma * v_0 * cos(particles[i].phi); // Update ax
                particles[i].ay += my_gamma * v_0 * sin(particles[i].phi); // Update ay
            }
        }
    }
    // Calculate and store the average potential energy per particle for a given timestep in the data set of the last particle in the timestep
    particles[numParticles - 1].PE = totalPotentialEnergy / particles.size();
}

// Function to update the positions and velocities using the Velocity Verlet method
void updatePositionsAndVelocities(vector<Particle>& particles, double L, double dt, vector<double>& vNc, vector<double>& vsum_term) {

    double k_B=1.0;
    double half_dt = 0.5 * dt;
    double var= sqrt(2.0 * my_gamma * k_B * temperature);
    // Generate random number generator for Gaussian white noise
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> noiseDist(0.0, 1.0);
    normal_distribution<double> Lambda_i(0.0, sqrt(2.0 * D_R));

    // Initialize vectors to store variables related to vision calculations and noise
    vector<vector<double>> noise(
        numParticles, vector<double>(
            4, 0.0
            )
        );

    int i=0;
    for (auto& particle : particles) {
        //Saves Eta and Zeta noise terms in each direction for each particle to be used in thsi and the next for loop
        double noise_eta = noiseDist(gen);
        double noise_zeta = noiseDist(gen);

        // Update velocities (half step)
        // Save the current particle's velocity and acceleration
        double dphi_i = 0.0;
        double v_x = particle.vx; double a_x = particle.ax;
        double v_y = particle.vy; double a_y = particle.ay;
        double v_i=sqrt(v_x*v_x+v_y*v_y);
        if (v_i==0) continue;
        noise[i][0] = noise_zeta*v_x/v_i;
        noise[i][1] = noise_zeta*v_y/v_i;
        noise[i][2] = noise_eta*v_x/v_i;
        noise[i][3] = noise_eta*v_y/v_i;

        // Update velocities using half step
        particle.vx += half_dt * (a_x - my_gamma * v_x) + 0.5 * sqrt(dt) * var * noise[i][0] - 0.125 * dt * dt * my_gamma * (a_x - my_gamma * v_x) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][0] + (1.0/sqrt(3.0)) * noise[i][2]);
        particle.vy += half_dt * (a_y - my_gamma * v_y) + 0.5 * sqrt(dt) * var * noise[i][1] - 0.125 * dt * dt * my_gamma * (a_y - my_gamma * v_y) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][1] + (1.0/sqrt(3.0)) * noise[i][3]);        
        // Update positions
        particle.x += particle.vx * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * noise[i][2];
        particle.y += particle.vy * dt + sqrt(dt) * dt * var * (0.5 / sqrt(3.0)) * noise[i][3];

        // Apply periodic boundary conditions
        particle.x -= floor(particle.x / L) * L;
        particle.y -= floor(particle.y / L) * L;
        i+=1;

    }
    // Update forces and calculate new potential energy
    LJandVissionCone(particles, vNc, vsum_term, L);

    i=0;
    for (auto& particle : particles) {
        // Save the current particle's new acceleration
        double dphi_i = 0.0;
        double v_xn = particle.vx; double a_xn = particle.ax;
        double v_yn = particle.vy; double a_yn = particle.ay;

        // Update velocities (ful step)
        particle.vx += half_dt * (a_xn - my_gamma * v_xn) + 0.5 * sqrt(dt) * var * noise[i][0] - 0.125 * dt * dt * my_gamma * (a_xn - my_gamma * v_xn) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][0] + (1.0/sqrt(3.0)) * noise[i][2]);
        particle.vy += half_dt * (a_yn - my_gamma * v_yn) + 0.5 * sqrt(dt) * var * noise[i][1] - 0.125 * dt * dt * my_gamma * (a_yn - my_gamma * v_yn) - 0.25 * sqrt(dt) * dt * my_gamma * var * (0.5 * noise[i][1] + (1.0/sqrt(3.0)) * noise[i][3]);
        i+=1;
    }
}

//Snippet to initialize particle possition in a cubic lattice
void initializeSystem(vector<Particle>& particles, double L, double p) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> uniform(-M_PI,M_PI);

    int particlesPerDimension = ceil(sqrt(numParticles));
    double spacing = L/double(particlesPerDimension);

    particles.resize(numParticles);

    // Calculate the number of particles with ID "1"
    int numParticles1 = round(numParticles * p);

    // Create a vector of particle IDs
    vector<int> particleIDs(numParticles, 2);
    fill(particleIDs.begin(), particleIDs.begin() + numParticles1, 1);

    // Shuffle the particle IDs
    shuffle(particleIDs.begin(), particleIDs.end(), gen);

    int index = 0;

    // Loop over each dimension of the cubic grid (x, y, z)
    for (int x = 0; x < particlesPerDimension; ++x) {
        for (int y = 0; y < particlesPerDimension; ++y) {
            // Check if all particles have been placed, and break the loop if so
            if (index >= numParticles) break;

            // Assign the particle ID and mass
            particles[index].particleID = particleIDs[index];
            particles[index].m = 1.0;
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

void saveParticleData( vector< vector< vector<double>>>& tensor, vector<Particle> particles, int step) {
    // Get the number of time steps, particles, and dimensions
    int numParticles = static_cast<int>(particles.size());

    int jsonStep=int(floor(double(step)/double(dataCompression)));
    // Fill the tensor with particle positions
    for (int p = 0; p < numParticles; ++p) {
        tensor[jsonStep][p][0] = particles[p].particleID; // Particle ID
        tensor[jsonStep][p][1] = particles[p].m; // Potential Energy
        tensor[jsonStep][p][2] = particles[p].x; // x-coordinate
        tensor[jsonStep][p][3] = particles[p].y; // y-coordinate
        tensor[jsonStep][p][4] = particles[p].vx; // x-velocity
        tensor[jsonStep][p][5] = particles[p].vy; // y-velocity
        tensor[jsonStep][p][6] = particles[p].phi; // Vission Cone direction
        if (p==numParticles-1) tensor[jsonStep][p][7] = particles[p].PE; // y-velocity
    }

}

void Notifier(int codeVersion, int currentVersion){
    //Plots the Energy.dat file using plotityy script and feeds it the used parameters which the script uses for naming the created html and png file
    string cd ="python Notifier.py "+ to_string(codeVersion) + " " + to_string(currentVersion);
    system(cd.c_str());
}

void createDirectory(const string& dir) {
    // program is being compiled on a Windows system. _mkdir function from the direct.h library to create the directory
    #ifdef _WIN32
        _mkdir(dir.c_str());
    // compiled on a Unix-based system. mkdir function from the sys/stat.h library to create the directory. The second argument 0777 sets the permissions of the directory 
    #else
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
    myfile << "theta,"<<theta*180/M_PI<<"\n";
    myfile << "R_0,"<<R_0<<"\n";
    myfile << "D_R,"<<D_R<<"\n";
    myfile << "D_T,"<<temperature/my_gamma<<"\n";
    myfile << "Pe,"<<(epsilon/temperature)-1.0<<"\n";
    myfile << "v_0,"<<v_0<<"\n";
    myfile << "Omega,"<<Omega<<"\n";
    myfile.close();
}
bool directoryExists(const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}
string dirlocater(){
    // Get the current time
    time_t now = time(0);
    char* dt = ctime(&now);

    // Convert ctime output to string and remove newline character
    string str_dt(dt);
    str_dt.erase(remove(str_dt.end()-1, str_dt.end(), '\n'), str_dt.end());

    // Replace spaces and colons with underscores
    replace(str_dt.begin(), str_dt.end(), ' ', '_');
    replace(str_dt.begin(), str_dt.end(), ':', '_');
    string folderName = "Runs/"+str_dt;

    if (!directoryExists(folderName)) {
        mkdir(folderName.c_str(), 0777);  // You can adjust the permissions as needed
        cout << "Directory created: " << folderName << endl;
    } 
    else {
        int count = 1;
        string numberedFolder = folderName + "_" + to_string(count);

        while (directoryExists(numberedFolder)) {
            count++;
            numberedFolder = folderName + "_" + to_string(count);
        }

        mkdir(numberedFolder.c_str(), 0777);  // Adjust permissions if needed
        folderName = numberedFolder;
        cout << "Directory created: " << numberedFolder << endl;
    }

    return folderName;
}

int main() {
    // // ask for a comment from the user
    // cout << "Please enter a comment for this simulation: ";
    // getline(cin, comment);
    // Set the number of threads for paralaization
    int numThreads = omp_get_max_threads(); // Set the desired number of threads
    omp_set_num_threads(numThreads);
    // get set number of threads from omp
    vector<Particle> particles(numParticles);
    vector<double> velocities;  // Empty velocity vector
    vector<double> vNc(numParticles, 0.0);
    vector<double> vsum_term(numParticles, 0.0);
    vector<vector<vector<double>>> tensor(
        int(floor(double(numSteps)/double(dataCompression))), vector<vector<double>>(
            numParticles, vector<double>(
                8, 0.0
            )
        )
    );
    //Calculating Packing fraction
    double packingFraction = (numParticles*M_PI*(sigma/2)*(sigma/2))/(boxSize*boxSize);
    cout << "Packing Fraction is: " <<packingFraction<<" at theta ="<<theta*180/M_PI <<endl;

    string location = dirlocater();
    save_param(location);
    
    // Main simulation loop
    initializeSystem(particles, boxSize, particleRatio); // Initialize the particle positions
    for (int step = 0; step < numSteps; ++step) {
        // Update positions and velocities
        updatePositionsAndVelocities(particles, boxSize, timestep, vNc, vsum_term);

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
    ofstream file(location+"/particle_positions.json");

    if (file.is_open()) {
        // Serialize the tensor to JSON
        nlohmann::json jsonData(tensor);

        // Write the JSON data to the file
        file << jsonData.dump(4); // The argument sets the indentation for readability
        file.close();
    } 
    else {cerr << "Unable to open the JSON file for writing." << endl;}

    // Print a success message
    cout << "Simulation completed successfully. Starting the Plotting and Notifying engines" << endl;


    // Notifier(codeVersion, currentVersion);

    return 0;
}
