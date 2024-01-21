#include<iostream>             
#include<string>              
#include<fstream>            
#include<cmath>                   
#include<cstdlib>             
#include<ctime>       
#include <vector>                                                                                           
using namespace std;

// Particle structure
struct Particle {
  double x, y, z, vx, vy, vz, ax1, ay1, az1, ax2, ay2, az2;
};

int numParticles = 864;
const double blx = 14.22757;
const double bly = 14.22757;
const double blz = 14.22757;
const double blxh = 0.5 * blx, blyh = 0.5 * bly, blzh = 0.5 * blz;
double temp = 1.0, dt = 0.005, sigma = 1, rcut = 2.5, epsilon = 1, mass = 1.00, gama = 0.01, k_B = 1.0;
double rcutsq = rcut * rcut;
int nstep = 20000, neq = 10000;

// Function prototypes

int main (int argc, char *argv[]);
void posvel_update(vector<Particle>& particles, double &ke,double &pe,double &tinst);
void calcforces(vector<Particle>& particles, double &pe) ;
void scales(vector<Particle>& particles,double &ke,double &tinst);
double gauss();




int main (int argc, char *argv[]){ 

  int iter;
  double ke,pe,kein,pein,etot,tinst;       
  srand(time(NULL));

  fprintf ( stderr, " --------- 1 ---------\n" ) ;
  fprintf( stderr, "collecting instantaneous datas\n" );

  //  Allocate memory.   
  vector<Particle> particles(numParticles);

  //Readng initial positions
  std::ifstream inputFile("fccdata.dat");
  if (!inputFile.is_open()){
    std::cerr << "Unable to open input file." << std::endl;
    return 1;
  }
  for (int i = 0; i < numParticles; ++i) {
    Particle particle;
    inputFile >> i >> particle.x >> particle.y >> particle.z;
    particles[i] = particle;
  }
  inputFile.close();

  
  // Save the values to a file
  std::ofstream outputFile("output.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Unable to open output file." << std::endl;
    return 1;
  }
  for (const Particle& particle : particles) {
    outputFile << "x: " << particle.x << " y: " << particle.y << " z: " << particle.z << std::endl;
  }
  outputFile.close();

  //Initializing velocity
  for(auto& particle : particles) {
  
    particle.vx = 0; // MB DIST...
    particle.vy = 0;
    particle.vz = 0;

    particle.ax2 = 0;
    particle.ay2 = 0;
    particle.az2 = 0;
    
  }  
  FILE *simulation;
  simulation=fopen("results1.dat","w"); 
                            
  for(iter=0;iter<=nstep;iter++) {
    posvel_update(particles,ke,pe,tinst); 
    if(iter<=neq) scales(particles,ke,tinst);

    kein = ke/double(particles.size());
    pein = pe/double(particles.size());
    etot = kein + pein;
    fprintf(simulation," %d\t%.8f\t%.8f\t%.8f\t%.8f\n",iter,pein,kein,etot,tinst) ;  
  } 


  fprintf ( stderr, " --------- 2 ---------\n" ) ;
  fprintf( stderr, "COMPLETED\n" );

  return 0;

}

void posvel_update(vector<Particle>& particles,double &ke,double &pe,double &tinst) { 


  double var = sqrt(2.0 * gama * k_B * temp);
  double c3 = 2/(3*float(particles.size()));

  
  double ctx,cty,ctz; 
  for (auto& particle : particles) {

    double noise_xi_x = gauss();
    double noise_eta_x= gauss();

    double noise_xi_y = gauss();
    double noise_eta_y= gauss();

    double noise_xi_z = gauss();
    double noise_eta_z= gauss(); 

    ctx = 0.5 * dt *dt * (particle.ax1 - gama * particle.vx) + var * sqrt(dt) * dt * (0.5 * noise_xi_x + (1/sqrt(3)) * noise_eta_x) ;
	  cty = 0.5 * dt *dt * (particle.ay1 - gama * particle.vy) + var * sqrt(dt) * dt * (0.5 * noise_xi_y + (1/sqrt(3)) * noise_eta_y) ;
	  ctz = 0.5 * dt *dt * (particle.az1 - gama * particle.vz) + var * sqrt(dt) * dt * (0.5 * noise_xi_z + (1/sqrt(3)) * noise_eta_z) ;

	  double diffx = particle.vx * dt + ctx;
	  double diffy = particle.vy * dt + cty;
	  double diffz = particle.vz * dt + ctz;

    particle.x = particle.x + diffx;
    particle.y = particle.y + diffy;
    particle.z = particle.z + diffz;    

    particle.ax2 = particle.ax1;
    particle.ay2 = particle.ay1;
    particle.az2 = particle.az1;
  }
  
  for (auto& particle : particles) {
    particle.x -= floor(particle.x / blx) * blx;
    particle.y -= floor(particle.y / bly) * bly;
    particle.z -= floor(particle.z / blz) * blz;
  }
  
  calcforces(particles,pe); 

  for (auto& particle : particles) {
    double noise_xi_x = gauss();
    double noise_eta_x= gauss();
    
    double noise_xi_y = gauss();
    double noise_eta_y= gauss();
    
    double noise_xi_z = gauss();
    double noise_eta_z= gauss();
        
    particle.vx = particle.vx + 0.5 * dt * (particle.ax1 + particle.ax2) - dt * gama * particle.vx + var * sqrt(dt) * noise_xi_x - gama * ctx;
    particle.vy = particle.vy + 0.5 * dt * (particle.ay1 + particle.ay2) - dt * gama * particle.vy + var * sqrt(dt) * noise_xi_y - gama * cty;
    particle.vz = particle.vz + 0.5 * dt * (particle.az1 + particle.az2) - dt * gama * particle.vz + var * sqrt(dt) * noise_xi_z - gama * ctz;
  } 

  ke = 0;  
      
  for(int i = 0;i < particles.size() ;i++)   ke = ke + (particles[i].vx * particles[i].vx  + particles[i].vy * particles[i].vy  + particles[i].vz * particles[i].vz);
  ke = ke*mass*0.5;
  tinst = ke*c3;
  return ;               
} 





/***********************************************************************************************/

void calcforces(vector<Particle>& particles, double &pe){
  for(auto& particle : particles) {
    particle.ax1 =0;
    particle.ay1 =0;
    particle.az1 =0;
  }
  pe=0; 

  for(int i = 0;i < particles.size();i++) { //no of neighbours for ith particle
    for(int j=i+1;j<particles.size();j++) {

      double dx = particles[i].x - particles[j].x;      
      double dy = particles[i].y - particles[j].y;
      double dz = particles[i].z - particles[j].z;
      if (dx>=blxh)       dx -= blx;
      else if (dx<-blxh)  dx += blx;
          
      if (dy>=blyh)       dy -= bly;
      else if (dy<-blyh)  dy += bly;    
      
      if (dz>=blzh)       dz -= blz;
      else if (dz<-blzh)  dz += blz;     
      
      double dsq = dx*dx+dy*dy+dz*dz;       
        
      if(dsq <= rcutsq) {
        double dr = sqrt(dsq);    
        
        double r_inv = sigma / dr;
        double r_inv6 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv ;
        double r_inv12 = r_inv6 * r_inv6;
        double force = -48 * (epsilon / dr) * (r_inv12 - 0.5 * r_inv6);
      // Add Gaussian white noise to forces
        particles[i].ax1 -= force * (dx / double(dr));
        particles[i].ay1 -= force * (dy / double(dr));
        particles[i].az1 -= force * (dz / double(dr));
        particles[j].ax1 += force * (dx / double(dr));
        particles[j].ay1 += force * (dy / double(dr));
        particles[j].az1 += force * (dz / double(dr));
        
        double potential = 4 * epsilon * (r_inv12 - r_inv6);
        pe += potential;         
      } //if loop
    }   //j loop
  }     //i loop
  return;
}



void scales(vector<Particle>& particles,double &ke,double &tinst){
  double ratio, c3 = 2/(3*float(particles.size()));
  ke=0;
  for(int i = 0;i < particles.size();i++)   ke = ke + (pow(particles[i].vx,2)+pow(particles[i].vy,2)+pow(particles[i].vz,2)); 
  ke=0.5*ke*mass;
  tinst=ke*c3;
  ratio=sqrt(temp/tinst);

  for(auto& particle : particles) {
    particle.vx = particle.vx * ratio;
    particle.vy = particle.vy * ratio;
    particle.vz = particle.vz * ratio;
  }   
  ke=0;
  for(int i = 0;i < particles.size();i++)   ke = ke + (pow(particles[i].vx,2)+pow(particles[i].vy,2)+pow(particles[i].vz,2)); 
  ke=0.5*ke*mass;
  tinst=ke*c3;
  return ;
}



double gauss(){
  static bool available = false;
  static double gset;
  double fac, rsq, v1, v2;
  if (!available) {
    do{
      v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
      v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    available = true;

    return v2*fac;
  } 
  else{
  available = false;
  return gset;
  }
}
