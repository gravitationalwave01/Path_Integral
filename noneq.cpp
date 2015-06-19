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

double simpson38(std::vector<std::pair<double,double> >* func) // implementing simpson's 3/8's rule to integrate func
{

  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly 
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin() + 3; it < func->end(); it+=3)
    ans += 3.0/8.0 * delta * ( (it-3)->second + 3 * (it-2)->second + 3 * (it-1)->second + it->second ); //simpson's 3/8 rule -- integrates three intervals at a time

  
  int rem = func->size() % 3; // if we have leftover intervals
  std::vector<std::pair<double,double> >::iterator it = func->end();
  if (rem != 0)
    if (rem == 1)
      ans += 0.5 * delta * ( (it-1)->second + (it-2)->second ); //trapezoid integration for the last interval
    if (rem == 2)
      ans += delta * ( (it-3)->second + 4 * (it-2)->second + (it-1)->second ) / 3; //simpsons rule for the last two intervals

  return ans;

}

int main()
{
  //initialize some constants
  int dim = 1;
  int P = 1;
  double time = 1000; // total time to go from lambda = 0 to lambda = 1
  double deltaT = 0.05;


  //open the output files:
  std::ofstream fpos("positions.dat"); //for positions
  std::ofstream fvel("velocities.dat"); //for velocities
  std::ofstream fene("energy.dat"); //for energies (kinetic and potential)
  std::ofstream fpot("fpot.dat");   // for dU/dlambda
  std::ofstream dist("dist.dat"); //stores the distribution of velocities (should be an MB distribution)
  std::ofstream work("work.dat"); //stores the computed values for the work

  //create the two potentials
  HarmonicPotential hp1(1.0,3.0);
  HarmonicPotential hp2(1.0,5.0);

  //Initialize the random number generators
  std::default_random_engine generator;
  double a = 1; // a^2 = kT/m
  std::normal_distribution<double> distribution(0,a);


  //Initialize a vector that will store an ensemble of ring polymers
  int numPolymers = 100000;
  std::vector<RP> ensemble;
  double tmp = 0;

  
  
  //initialize the ensemble with RP positions and velocities randomly from this distribution
  std::cout << " **** STARTING RP INIT  ****" << std::endl;
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
	velocities.push_back(distribution(generator));
	positions.push_back(0);
      }	
    }
    
    //create a polymer with these initial conditions and store in the ensemble 
    ensemble.push_back(RP(dim,P,&positions,&velocities,&masses));
  }
  
  
  //we have the ensemble initialized -- now propogate the RP in time so that they equilibrate.
  std::vector<double> curForces; 
  std::vector<double> oldForces; 

  double kinetic, potential, pot;//pot stores dU/dlambda
  bool first = true;

  

  //as we vary lambda, the potential changes and we need to wait for the system to equlibrate. 
  //after equlibrating we can sample the potential energy
  int numSamples = 0;
  double timeSinceLast = 0; //time since the last collision with the heat bath
  int numCols = 0; //number of collisions
  double prefactor = 0.5 * hp1.mass * (hp2.omega*hp2.omega - hp1.omega*hp1.omega);
  // Smoothly vary the potential-flipping parameter:
  double dlambda = deltaT/time;
  double exp_work;
  double progress = 0;

  pot = 0;
  numSamples = 0;
  int counter = 0;
  double lambda;
  for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
  {
    std::vector<std::pair<double, double> > integrand; //will store the values (lambda, dV/dlambda )
    lambda = 0;
    for (double curTime = 0; curTime < time; curTime+=deltaT)
    {
      //update forces:
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
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	  curForces[curBead*dim+j] = (1-lambda)*hp1.getForce(myRP->beads[curBead].position[j]) + lambda*hp2.getForce(myRP->beads[curBead].position[j]);
      
      //update velocities
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	  myRP->beads[curBead].velocity[j] += deltaT * 0.5 / myRP->beads[curBead].mass * (oldForces[curBead*dim+j]+curForces[curBead*dim+j]);
      
      //compute current potential and kinetic energies
      for (int curBead = 0; curBead < myRP->P; curBead++)
	for (int j=0;j<dim;j++)
	{
	  kinetic += 0.5 * myRP->beads[curBead].mass * myRP->beads[curBead].velocity[j] * myRP->beads[curBead].velocity[j];
	  potential += ( (1-lambda)*hp1.getPotential(myRP->beads[curBead].position[j]) + lambda*hp2.getPotential(myRP->beads[curBead].position[j]) );
	}


      // sample dU/dlambda with some probability.
      // to sample dU/dlambda we really sample <x>, plug this x into (U_{lambda+dlambda/2) - U_{lambda-dlambda}) / 2
      //      if (distribution(generator) >= 0 && timeSinceLast > 0) // 
      pot = 0;
      for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	for (int j=0;j<dim;j++)
	{
	  pot += hp2.getPotential(it->position[j]) - hp1.getPotential(it->position[j]);
	  //	  dist << lambda << " " << it->velocity[j] << " " << it->position[j] << std::endl;
	}
      integrand.push_back(std::pair<double,double>(curTime,pot)); 
      lambda += dlambda;

    } //end of time loop

    counter++;
    tmp = 1/time * simpson38(&integrand);
    if ((double)counter/numPolymers > progress)
    {
      std::cout << progress*100 << "% done " << std::endl;
      progress += 0.1;
    }
    work << tmp << std::endl;
    exp_work += exp(1/a*1/a * tmp);//assumes that mass=1 therefore a^2 = T
  } // end of ensemble loop

  //close all output files:
  fpos.close();
  fvel.close();
  fene.close();
  fpot.close();


  //output the obtained values of dU/dlambda
  //  std::cout << "Here are the values of the potential" << std::endl;
  //  for (std::vector<std::pair<double,double> >::iterator it = integrand.begin(); it < integrand.end(); it++)
  //    std::cout << it->first << " " << it->second << std::endl;
  //  std::cout << "Done outputting the values of the potential" << std::endl;

  //  std::cout << " NUMBER OF COLLISIONS: " << numCols << std::endl;
  //  std::vector<std::pair<double,double> > test;
  std::cout << "FREE ENERGY DIFFERENCE IS : " << std::endl;
  std::cout << a*a*log(exp_work/numPolymers) << std::endl;
  
};
