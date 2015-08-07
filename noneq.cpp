// Now I am trying to implement a thermalized (isokinetic) dynamics.


#include <utility>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>

//#define pb push_back;

unsigned long long rdtsc()
{
  unsigned int lo,hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return ((unsigned long long)hi << 32) | lo;
}

class QuarticPotential
{
 public:
  double k;

  QuarticPotential(double k_tmp)
  {
    k = k_tmp;
  }

  double getPotential(double pos)
  {
    return 0.5 * k * pos * pos * pos * pos;
  }

  double getForce(double pos)
  {
    return - 2 * k * pos * pos * pos;
  }
};


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

  RP(int d, int numBeads, std::vector<double>* positions, std::vector<double>* velocities, std::vector<double>* masses)
  {
    P = numBeads;
    dim = d;
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


int main()
{
  //initialize some constants
  int dim = 1;
  int P = 1;
  double timevec[5] = {1.0,3.0,10.0,30.0,100.0}; // total time to go from lambda = 0 to lambda = 1
  double time;
  double deltaT = 0.01;
  double collFreqVec[5] = {0.05,0.01,0.005,0.001,0.0005}; //strength of interaction with bath
  double collFreq;

  //  double a = 1.2247; // a^2 = kT/m
  double ibeta = 1.5; //inverse beta, ibeta = kT

  //create the two potentials
  double omega1 = 1;
  double omega2 = 2;
  HarmonicPotential hp1(1.0,omega1);
  HarmonicPotential hp2(1.0,omega2);
  double answer = ibeta * log(omega2/omega1);
  std::ofstream tmp_out("tmp.dat");

  int numTrials = 25;
  std::vector<double> deltaF;
  for (int trial = 0; trial < numTrials; trial++)
  {
    std::cout << "starting trial " << trial << " of " << numTrials << std::endl;
    time = timevec[trial%5];
    collFreq = collFreqVec[trial/5];
    tmp_out << time << " " << collFreq << " " ;
  //open the output files:
  //  std::ofstream fpos("positions.dat"); //for positions
  //  std::ofstream fvel("velocities.dat"); //for velocities
  std::ofstream fene("energy.dat"); //for energies (kinetic and potential)
  //  std::ofstream fpot("fpot.dat");   // for dU/dlambda
  std::ofstream dist("dist.dat"); //stores the distribution of velocities (should be an MB distribution)
  std::ofstream work("work.dat"); //stores the computed values for the work

  //Initialize the random number generators
  std::default_random_engine generator(rdtsc());
  std::normal_distribution<double> pdist(0,sqrt(ibeta)); //momentum distribution
  std::normal_distribution<double> xdist(0,1/omega1*sqrt(ibeta)); //position distribution
  std::uniform_real_distribution<double> unif(0.0,1.0);


  //Initialize a vector that will store an ensemble of ring polymers
  int numPolymers = 10000;
  std::vector<RP> ensemble;
  double tmp = 0;


  //initialize the ensemble with RP positions and velocities randomly from this distribution
  for (int i = 0; i < numPolymers; i++)
  {
    std::vector<double> positions;
    std::vector<double> velocities;
    std::vector<double> masses;
    
    for (int curBead = 0; curBead < P; curBead++)
    {
      masses.push_back(1);
      for (int curDim = 0; curDim < dim; curDim++)
      {
	velocities.push_back(pdist(generator));
	positions.push_back(xdist(generator));
      }	
    }
    
    //create a polymer with these initial conditions and store in the ensemble 
     ensemble.push_back(RP(dim,P,&positions,&velocities,&masses));
  }
  
  
  //we have the ensemble initialized
  std::vector<double> curForces; 
  std::vector<double> oldForces; 

  double pot;//pot stores dU/dlambda
  bool first = true;

  

  double dlambda = deltaT/time;
  double exp_work=0;
  double progress = 0;
  double progress2 = collFreq;
  unsigned long int numcoll = 0;


  pot = 0;
  int counter = 0;
  double lambda =0;
  //  std::vector<std::pair<double,double> > fart;
  //  std::vector<std::vector<std::pair<double, double> > > integrand(numPolymers, fart); 
  std::vector<std::vector<std::pair<double, double> > > integrand; 
  for (int i=0;i<numPolymers;i++)
  {
    std::vector<std::pair<double,double> > blah;
    integrand.push_back(blah);
  }  

  //the size of integrand is numPolymers
  std::vector<std::vector<std::pair<double, double> > >::iterator it_integrand; 

  for (double curTime = 0; curTime < time; curTime+=deltaT)
  {
    
    //add in an Andersen thermostat

    //    if ((double)lambda > progress2)
    {
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	  if (unif(generator) < deltaT * collFreq) 
	    for (int j=0;j<dim;j++)
	      it->velocity[j] = pdist(generator);
    }
    
    //save the current total energy of the first RP to a file
    /*
    for (int i = 0; i < ensemble[0].P; i++)
      for (int j=0; j < dim; j++)
	fene << curTime << " " << 0.5*ensemble[0].beads[i].velocity[j]*ensemble[0].beads[i].velocity[j] + 
	  lambda*hp2.getPotential(ensemble[0].beads[i].position[j]) + (1-lambda) * hp1.getPotential(ensemble[0].beads[i].position[j]) << std::endl;
    */

    it_integrand = integrand.begin();
    for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
    {
      //update forces:
      //for RP: need to add in inter-bead forces
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	{
	  if (first)
	  {
	    tmp = (1-lambda)*hp1.getForce(myRP->beads[curBead].position[j]) + lambda*hp2.getForce(myRP->beads[curBead].position[j]);
	    curForces.push_back(tmp);
	    oldForces.push_back(tmp);
	    first = false;
	  } else
	  {
	    tmp = (1-lambda)*hp1.getForce(myRP->beads[curBead].position[j]) + lambda*hp2.getForce(myRP->beads[curBead].position[j]);
	    curForces[curBead*dim+j] = tmp;
	    oldForces[curBead*dim+j] = tmp;
	  }
	}
	
      //update positions
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	  myRP->beads[curBead].position[j] += deltaT * myRP->beads[curBead].mass * myRP->beads[curBead].velocity[j] 
	    + deltaT * deltaT * 0.5 / myRP->beads[curBead].mass * curForces[curBead*dim+j];
      
      //update current forces while keeping the old ones
      //TODO: add in inter-bead forces
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	  curForces[curBead*dim+j] = (1-lambda)*hp1.getForce(myRP->beads[curBead].position[j]) + lambda*hp2.getForce(myRP->beads[curBead].position[j]);
      
      //update velocities
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	  myRP->beads[curBead].velocity[j] += deltaT * 0.5 / myRP->beads[curBead].mass * (oldForces[curBead*dim+j]+curForces[curBead*dim+j]);
      
      // sample dU/dlambda with some probability.
      // to sample dU/dlambda we really sample x, plug this x into (U_{lambda+dlambda/2) - U_{lambda-dlambda}) / 2
      pot = 0;
      for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	for (int j=0;j<dim;j++)
	    pot += ( hp2.getPotential(it->position[j]) - hp1.getPotential(it->position[j]) ) ;

      it_integrand->push_back(std::pair<double,double>(curTime,pot)); 
      it_integrand++;

    } //end of ensemble loop


    lambda += dlambda;
    counter++;
  } // end of time loop

  //close all output files:
  //  fpos.close();
  //  fvel.close();
  fene.close();
  //  fpot.close();
  double tmp_avg=0;
  for (it_integrand = integrand.begin();it_integrand<integrand.end();it_integrand++)
  {
    tmp = 1/time * simpson38(&*it_integrand);
    work << tmp << std::endl;
    exp_work += exp(-1/ibeta * tmp); //assumes that mass=1 therefore a^2 = kT
    tmp_avg += tmp;
  }

  //  std::cout << "W^a is " << tmp_avg / numPolymers << std::endl;
  //  std::cout << "W^z is " << -ibeta*log(exp_work/numPolymers) << std::endl;
  tmp_out << tmp_avg / numPolymers << " " << -ibeta*log(exp_work/numPolymers) << std::endl;

  deltaF.push_back(-ibeta*log(exp_work/numPolymers));
  //  std::cout << "Ratio of collisions with the bath : "<< (1 / collFreq) / (time/deltaT) << std::endl;
  }//end of trial loop

  double finalanswer = 0;
  for (int i=0;i<numTrials;i++)
    finalanswer += deltaF[i];

  finalanswer/=numTrials;

  double stdDev = 0;
  for (int i=0;i<numTrials;i++)
    stdDev += (deltaF[i]-finalanswer)*(deltaF[i]-finalanswer);

  stdDev /=numTrials;
  stdDev = sqrt(stdDev);

  //compute and print the free energy difference
  /*
  std::cout << "COMPUTED FREE ENERGY DIFFERENCE via " << numTrials << " trials IS : " << finalanswer << std::endl;
  std::cout << "EXPECTED " << answer << std::endl;
  std::cout << "ERROR " << ((finalanswer-answer)/answer)*100 << "%" << std::endl;
  std::cout << "STD DEV of delta F values is " << stdDev << std::endl;
  */
};
