//This code computes equilibrium averages of quantum potentials using path integrals


#include <utility>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>
#include <sstream>
//#define pb push_back;

unsigned long long rdtsc()
{
  unsigned int lo,hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return ((unsigned long long)hi << 32) | lo;
}

class HarmonicPotential
{
 public:
  double mass;
  double omega;

  HarmonicPotential(double m, double o)
  {
    mass = m;
    omega = o;
  }

  double getPotential(double pos)
  {
    return 0.5 * mass * omega * omega * pos * pos;
  }

  double getForce(double pos)
  {
    return -mass * omega * omega * pos;
  }
};


class Bead
{
 public:
  int dim; //dimension of the system
  double* position;
  double* velocity;
  double mass;
  
  Bead(std::vector<double>::iterator pos, int d, std::vector<double>::iterator vel, double m)
  {
    dim = d;
    position = new double[dim];
    velocity = new double[dim];
    for (int i=0;i<dim;i++)
    {
      position[i] = pos[i];
      velocity[i] = vel[i];
    }
    mass = m;
  }


  void printMe()
  {
    std::cout << "bead position: ";
    for (int i=0;i<dim;i++)
      std::cout << position[i] << " ";
    std::cout << std::endl;
    
    std::cout << "bead velocity: ";
    for (int i=0;i<dim;i++)
      std::cout << velocity[i] << " ";
    std::cout << std::endl;
  }
};


class RP
{
 public:
  int dim; //dimension of the system
  std::vector<Bead> beads;
  int P; //number of beads
  double omega;

  RP(int d, int numBeads, std::vector<double>* positions, std::vector<double>* velocities, std::vector<double>* masses, double o)
  {
    P = numBeads;
    dim = d;
    omega = o;

    for (int i =0; i < numBeads; i++)
    {
      beads.push_back(Bead(positions->begin()+i*dim, dim, velocities->begin()+i*dim, masses->at(i)));
    }
  }

  void printMe()
  {
    std::cout << "DETAILS OF THE RING POLYMER::" << std::endl;
    std::cout << "POSITIONS:     VELOCITIES:" << std::endl;
    for (int i = 0; i < P; i++)
      for (int j=0;j<dim;j++)
      std::cout << beads[i].position[j] << "                       " << beads[i].velocity[j] << std::endl;
  }
};


//this function integrates *func using the traezoidal method.
double integrate(std::vector<std::pair<double,double> >* func)
{
  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly and constantly
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin(); it+1 < func->end(); it++)
    ans += delta * ( (it+1)->second + it->second) / 2;
  
  return ans;
}

double simpson(std::vector<std::pair<double,double> >* func) // implementing simpson's rule to integrate func
{
  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly 
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin() + 2; it < func->end(); it+=2)
  {
    ans += delta * ( (it-2)->second + 4 * (it-1)->second + it->second ) / 3; //simpson's 3/8 rule -- integrates three intervals at a time
    //    std::cout << ans << std::endl;
  }

  int rem = func->size() % 2; // if we have leftover intervals
  std::vector<std::pair<double,double> >::iterator it = func->end();
  //  std::cout << "rem is " << rem << " " << func->size() << std::endl;
  if (rem == 0)
    ans += 0.5 * delta * ( (it-1)->second + (it-2)->second ); //trapezoid integration for the last interval

  return ans;
}

double simpson38(std::vector<std::pair<double,double> >* func) // implementing simpson's 3/8's rule to integrate func
{

  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly 
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin() + 3; it < func->end(); it+=3)
  {
    ans += 3.0/8.0 * delta * ( (it-3)->second + 3 * (it-2)->second + 3 * (it-1)->second + it->second ); //simpson's 3/8 rule -- integrates three intervals at a time
    //    std::cout << ans << std::endl;
  }
  
  int rem = func->size() % 3; // if we have leftover intervals
  std::vector<std::pair<double,double> >::iterator it = func->end();
  //  std::cout << "rem is " << rem << " " << func->size() << std::endl;
  if (rem != 1)
    if (rem == 2)
      ans += 0.5 * delta * ( (it-1)->second + (it-2)->second ); //trapezoid integration for the last interval
    else if (rem == 0)
      ans += delta * ( (it-3)->second + 4 * (it-2)->second + (it-1)->second ) / 3; //simpson's rule for last two intervals

  return ans;
}


//propogates a single ring polymer in time by a single time step deltaT
double velocityVerlet(RP* myRP, double deltaT, HarmonicPotential* hp)
{
  std::vector<double> curForces;
  std::vector<double> oldForces;
  curForces.resize(myRP->P);
  oldForces.resize(myRP->P);
  double dim = myRP->dim;
  double tmp=0;
  double omega = myRP->omega;
  int P = myRP->P;
  

  //update forces:
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
    {
      tmp =  omega*omega*myRP->beads[curBead].mass * 
	( (myRP->beads[((curBead-1)%P+P)%P].position[j] - myRP->beads[curBead%P].position[j] )  //shitty C++ modular arithmetic: (-1%4)=-1, not 3 :( 
	  + (myRP->beads[(curBead+1)%P].position[j] - myRP->beads[curBead%P].position[j] ) );
      
      tmp += hp->getForce(myRP->beads[curBead].position[j]) / myRP->P; //accounts for external potential

      curForces[curBead*dim+j] = tmp;
      oldForces[curBead*dim+j] = tmp;
    }  

  //update positions
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
      myRP->beads[curBead].position[j] += deltaT * myRP->beads[curBead].mass * myRP->beads[curBead].velocity[j] 
	+ deltaT * deltaT * 0.5 / myRP->beads[curBead].mass * curForces[curBead*dim+j];
  
  //update current forces while keeping the old ones
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
    {
      tmp =  omega*omega*myRP->beads[curBead].mass * 
	( (myRP->beads[((curBead-1)%P+P)%P].position[j] - myRP->beads[curBead%P].position[j] ) 
	  + (myRP->beads[(curBead+1)%P].position[j] - myRP->beads[curBead%P].position[j] ) );
      
      tmp += hp->getForce(myRP->beads[curBead].position[j]) / myRP->P; //accounts for external potential

      curForces[curBead*dim+j] = tmp;
    }
  
  //update velocities
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
      myRP->beads[curBead].velocity[j] += deltaT * 0.5 / myRP->beads[curBead].mass * (oldForces[curBead*dim+j]+curForces[curBead*dim+j]);

  /*
  double kinetic = 0;
  double extpot = 0;
  double beadpot = 0; 
  //compute the energy of the system
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
    {
      kinetic += 0.5 * myRP->beads[curBead].velocity[j] * myRP->beads[curBead].velocity[j]; // mv^2/2
      extpot += hp->getPotential(myRP->beads[curBead].position[j])/myRP->P; // external potential
      //      std::cout << "computing difference between " << curBead << " and " << (curBead+1) % myRP->P << std::endl;
      beadpot += 0.5 * omega * omega * pow(myRP->beads[curBead].position[j] - myRP->beads[(curBead+1) % myRP->P].position[j],2); //spring potential
    }
  
  fene2 << time << " " << kinetic << " " << extpot << " " << beadpot << std::endl;
  */
  return 0;
};





int main()
{
  //initialize some constants
  int dim = 1;
  int Pvec[4] = {4,8,16,32}; //number of blocks
  int P;
  double ibetavec[1] = {0.5}; // ibeta = 1/beta = kT
  double ibeta;
  double timevec[1] = {50}; // total time to go from lambda = 0 to lambda = 1
  double time;
  double deltaT = 0.01;
  double collFreq = 0.002; //frequency of collision with bath
  HarmonicPotential hp(1,1);


  //Initialize the random number generators
  std::default_random_engine generator(rdtsc());
  std::uniform_real_distribution<double> unif(0.0,1.0);

  int numTrials = 1;

  //open the output files:
  std::ofstream fpos("positions.dat"); //for positions
  //  std::ofstream fvel("velocities.dat"); //for velocities
  //  std::ofstream fene("energy.dat"); //for energies (kinetic and potential)
  std::ofstream fpot("fpot.dat");   // for dU/dlambda
  std::ofstream dist("dist.dat"); //stores the distribution of velocities (should be an MB distribution)
  std::ofstream work("work.dat"); //stores the computed values for the work
  int curoutinc = 4;
  int curout = curoutinc;
  int obs_trial = 0; //the trial about which we record dyanmics, work, distributions, etc.
  int collCount = 0; // counts the number of collisions with the heat bath

  for (int trial = 0; trial < numTrials; trial++)
  {
    collCount = 0;
    time = timevec[0];
    std::cout << "starting trial " << trial+1 << " of " << numTrials << std::endl;

    P = Pvec[trial];
    double ibeta = ibetavec[0];
    double Pibeta = P * ibeta; // beta_P = beta/P => ibeta_P = P*ibeta 
    double Pomega = Pibeta; // spring stiffness
    double tmp = 0;
    std::normal_distribution<double> pdist(0,sqrt(Pibeta)); //boltzmann distribution for momenta
    tmp = sqrt(0.5 + 1/(exp(1/ibeta)-1));
    std::normal_distribution<double> xdist(0,tmp); //boltzmann distribution for positions



    //Initialize a vector that will store an ensemble of ring polymers
    int numPolymers = 500;
    std::vector<RP> ensemble;



    //initialize the ensemble with RP velocities randomly from this distribution
    //this equilibrates the initial state
    double comp;
    for (int i = 0; i < numPolymers; i++)
    {
      std::vector<double> positions;
      std::vector<double> velocities;
      std::vector<double> masses;

      comp = 0;
      for (int curBead = 0; curBead < P; curBead++)
      {
	masses.push_back(1.0);
	for (int curDim = 0; curDim < dim; curDim++)
	{
	  tmp = pdist(generator);
	  comp += tmp; //ignores dimensions > 1 !!!!!
	  velocities.push_back(tmp);
	  positions.push_back(xdist(generator)); //equilibrium position
	}	
      }
    
      for (int curBead = 0; curBead < P; curBead++)
	positions[curBead] -= comp / P;

      //create a polymer with these initial conditions and store in the ensemble 
      ensemble.push_back(RP(dim,P,&positions,&velocities,&masses,Pomega));
    }
  
    //we have the ensemble initialized with the proper kinetic energy, but we still need to get the right configuration
    //for this, we propogate in time for a bit:
    for (double time_t = 0; time_t < curoutinc*4 + 0.5; time_t += deltaT) //we don't care too much about energy conservation -- just get a distribution!
    {
      /*
      fpos << time_t << " " << ensemble[0].beads[0].position[0] 
	   << " " << ensemble[0].beads[1].position[0] 
	   << " " << ensemble[0].beads[2].position[0] 
	   << " " << ensemble[0].beads[3].position[0] << std::endl;
      */
      if ( time_t > curout ) 
      {
	std::stringstream ss;
	ss << "xdist" << curout/curoutinc << ".dat";
	std::ofstream xdist(ss.str().c_str());
	ss.str("");
	ss << "vdist" << curout/curoutinc << ".dat";
	std::ofstream vdist(ss.str().c_str());

	for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	  for (int curBead = 0; curBead < myRP->P; curBead++)
	  {
	    xdist << myRP->beads[curBead].position[0] << std::endl;
	    vdist << myRP->beads[curBead].velocity[0] << std::endl;
	  }

	curout+=curoutinc;
      }
      
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	velocityVerlet(&*myRP, deltaT, &hp);


      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	  if (unif(generator) < deltaT * collFreq) 
	  {
	    for (int j=0;j<dim;j++)
	      it->velocity[j] = pdist(generator); //only while the mass = 1
	    collCount++;
	  }
    }
    
    for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
      if (trial == obs_trial)
	dist << myRP->beads[0].position[0] << std::endl;
      
    //now the kinetic energies are all messed up, so re-assign them:
    for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
      for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	for (int j=0;j<dim;j++)
	  it->velocity[j] = pdist(generator); //only while the mass = 1
    
    
    double sample; //sample stores the quantity that I am sampling
    unsigned long numSamples = 0;
    double dlambda = deltaT/time;
    double exp_work=0;
    double lambda = 0;
    

    std::cout << " ..... INITIALIZED " << std::endl;
    //    std::cout << "time of simulation is " << time << ", deltaT is " << deltaT << std::endl;
    //    std::cout << "the number of collisions was " << collCount << ", which is once every " << ((time/deltaT)/collCount) << "time steps " << std::endl;

    lambda = 0;
    collCount = 0;
    for (double curTime = 0; curTime <= time; curTime+=deltaT)
    {
      
      if (trial == obs_trial)
	fpos << lambda << " " << ensemble[0].beads[0].position[0] 
	     << " " << ensemble[0].beads[1].position[0] 
	     << " " << ensemble[0].beads[2].position[0] 
      	     << " " << ensemble[0].beads[3].position[0] << std::endl;
      
      //Andersen thermostat
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	  if (unif(generator) < deltaT * collFreq) 
	  {
	    for (int j=0;j<dim;j++)
	      it->velocity[j] = pdist(generator); //only while the mass = 1
	    collCount++;
	    //	    std::cout << "collision!" << std::endl;
	  }


      //propagate in time
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	velocityVerlet(&*myRP,deltaT,&hp);

      // sample with some probability.
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	if (unif(generator) < 0.25)
	  for(int i = 0; i < myRP->P; i++)
	    for (int j = 0; j < dim; j++)
	    {
	      sample += pow(myRP->beads[i].position[j],2);
	      numSamples++;
	    }
	  
      
      lambda += dlambda;
    } // end of time loop
    

    std::cout << "time of simulation is " << time << ", deltaT is " << deltaT << std::endl;
    std::cout << "the total number of collisions was " << collCount << ", which is once every " << ((time/deltaT)/(collCount/P)) << "time steps PER BEAD" << std::endl;
    std::cout << "the expectation value is " << sample/numSamples << std::endl;

    //close all output files:
    //  fpos.close();
    //  fvel.close();
    //  fene.close();
    //  fpot.close();

  }//end of trial loop
  

};
