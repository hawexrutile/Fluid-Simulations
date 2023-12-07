
import math

epsilon = 1.0;           
sigma = 1.0;                 
numParticles = 625;                       
temperature = 1.0;              
tau=math.sqrt(sigma*sigma/temperature)
gamma=100/tau
DR=0.08/tau
DT= temperature/gamma
Pe=(epsilon/temperature)-1
v=Pe * DT / sigma
packingFraction = 0.6
# boxSize = math.sqrt((numParticles*math.pi*(sigma/2)*(sigma/2))/packingFraction);
boxSize = 100;
packingFraction = (numParticles*math.pi*(sigma/2)*(sigma/2))/(boxSize**2);
#print all the variables
print ("epsilon = ", epsilon)
print ("sigma = ", sigma)
print ("numParticles = ", numParticles)
print ("boxSize = ", boxSize)
print ("temperature = ", temperature)
print ("gamma = ", gamma)
print ("tau = ", tau)
print ("DR = ", DR)
print ("DT = ", DT)
print ("Pe = ", Pe)
print ("v = ", v)
print ("packingFraction = ", packingFraction)
print ("boxSize = ", boxSize)

