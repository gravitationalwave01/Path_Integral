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
    ans += delta * ( (it-2)->second + 4 * (it-1)->second + it->second ) / 3; //simpson's 3/8 rule -- integrates three intervals at a time


  int rem = func->size() % 2; // if we have leftover intervals
  std::vector<std::pair<double,double> >::iterator it = func->end();
  if (rem == 0)
    ans += 0.5 * delta * ( (it-1)->second + (it-2)->second ); //trapezoid integration for the last interval

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
  if (rem != 1)
    if (rem == 2)
      ans += 0.5 * delta * ( (it-1)->second + (it-2)->second ); //trapezoid integration for the last interval
    else if (rem == 0)
      ans += delta * ( (it-3)->second + 4 * (it-2)->second + (it-1)->second ) / 3; //simpson's rule for last two intervals

  return ans;
}


int main()
{
  int pspace_ctr = 0;
  std::vector<double> v_time{100000};
  std::vector<double> v_dlambda{0.02};
  std::vector<double> v_deltaT{0.01};
  std::ofstream result("./statistics/deltaF.dat");

  result << "The output is formatted as follows: " << std::endl;
  result << "pindex   time    deltaT    dlambda   trial  simp38   simp   trap " << std::endl;

  //Initialize the random number generator
  std::default_random_engine generator;
  double a = 1; // a^2 = kT/m
  std::normal_distribution<double> distribution(0,a);


  for (std::vector<double>::iterator vt_it = v_time.begin(); vt_it < v_time.end(); vt_it++)
  {
  for (std::vector<double>::iterator vl_it = v_dlambda.begin(); vl_it < v_dlambda.end(); vl_it++)
  {
  for (std::vector<double>::iterator vdt_it = v_deltaT.begin(); vdt_it < v_deltaT.end(); vdt_it++)
  {

  for (int trial = 0; trial < 5; trial++)
  {
    //    generator.seed(rdtsc());
  double time = *vt_it;
  double deltaT = *vdt_it;
  double dlambda = *vl_it;

  //initialize some constants
  int dim = 1;
  int P = 1;
  //  double time = 50000;
  //  double deltaT = 0.01;


  //open the output files:
  //  std::ofstream fpos("positions.dat"); //for positions
  //  std::ofstream fvel("velocities.dat"); //for velocities
  //  std::ofstream fene("energy.dat"); //for energies (kinetic and potential)

  std::stringstream fname("");
  fname  << "./statistics/fpot" << pspace_ctr << "_t" << trial << ".dat"; //name of output file.
  std::ofstream fpot(fname.str().c_str());   // file for storing dU/dlambda values
  //  std::ofstream dist("dist.dat"); //stores the distribution of velocities (should be an MB distribution)

  //create the two potentials
  HarmonicPotential hp1(1.0,3.0);
  HarmonicPotential hp2(1.0,12.0);
  double answer = log(12.0/3.0);


  //Initialize a vector that will store an ensemble of ring polymers
  int numPolymers = 1;
  std::vector<RP> ensemble;
  double tmp = 0;


  //initialize the ensemble with RP positions and velocities randomly from this distribution
  //  std::cout << " **** STARTING RP INIT  ****" << std::endl;
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
  std::vector<std::pair<double, double> > integrand; //will store the values (lambda, potential energy)
  double kinetic, potential, pot;//pot stores dU/dlambda
  bool first = true;

  

  //as we vary lambda, the potential changes and we need to wait for the system to equlibrate. 
  //after equlibrating we can sample the potential energy
  int numSamples = 0;
  double timeSinceLast = 0; //time since the last collision with the heat bath
  int numCols = 0; //number of collisions
  double prefactor = 0.5 * hp1.mass * (hp2.omega*hp2.omega - hp1.omega*hp1.omega);
  // Smoothly vary the potential-flipping parameter:
  //  double dlambda = 0.02;
  for (double lambda = 0; lambda <= 1.0001; lambda += dlambda)
  {
    pot = 0;
    numSamples = 0;
    //    std::cout << " ************** CURRENT LAMBDA ************* " << lambda << std::endl;    
    timeSinceLast = 100000; //to ensure that the first thing to happen is a collision with bath

    for (double curTime = 0; curTime < time; curTime+=deltaT)
    {
      //collide with the bath:
      if (timeSinceLast > 10)
      {
	for (std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	  for (int curBead = 0; curBead < myRP->P; curBead++)
	    for (int curDim = 0; curDim < myRP->dim; curDim++)
	    {
	      myRP->beads[curBead].velocity[curDim] = distribution(generator);
	    }
	timeSinceLast = 0;
	numCols++;
      }

      //      kinetic=0;
      //      potential=0;
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
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
	//	for (int curBead = 0; curBead < myRP->P; curBead++)
	//	  for (int j=0;j<dim;j++)
	//	  {
	//	    kinetic += 0.5 * myRP->beads[curBead].mass * myRP->beads[curBead].velocity[j] * myRP->beads[curBead].velocity[j];
	//	    potential += ( (1-lambda)*hp1.getPotential(myRP->beads[curBead].position[j]) + lambda*hp2.getPotential(myRP->beads[curBead].position[j]) );
	//	  }


      } //end of VV loop
      

      // sample dU/dlambda with some probability.
      // to sample dU/dlambda we really sample <x>, plug this x into (U_{lambda+dlambda/2) - U_{lambda-dlambda}) / 2
      if (distribution(generator) >= 0 && timeSinceLast > 5) // 
	for (std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	  for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	    for (int j=0;j<dim;j++)
	    {
	      pot += hp2.getPotential(it->position[j]) - hp1.getPotential(it->position[j]);
	      numSamples++;
	      //	      dist << lambda << " " << it->velocity[j] << " " << it->position[j] << std::endl;
	    }

      timeSinceLast += deltaT; 
   } //end of time loop

    integrand.push_back(std::pair<double,double>(lambda,pot/numSamples));
    //fpot << lambda << " " << pot/numSamples << " " << numSamples << std::endl;
  } // end of lambda loop

  //close all output files:
  //  fpos.close();
  //  fvel.close();
  //  fene.close();
  //fpot.close();
  //

  //output the obtained values of dU/dlambda
  //  std::cout << "Here are the values of the integrand:" << std::endl;
  for (std::vector<std::pair<double,double> >::iterator it = integrand.begin(); it < integrand.end(); it++)
    fpot << it->first << " " << it->second << std::endl;
  fpot.close();


  result << pspace_ctr << " " << time << " " << deltaT << " " << dlambda << " " << trial << " " << 
    std::abs(simpson38(&integrand)-answer)/answer*100 << " " << std::abs(simpson(&integrand)-answer)/answer*100 << " " << std::abs(integrate(&integrand)-answer)/answer*100 << std::endl;

  }
  //  result << std::endl;
  pspace_ctr++;
  }
  }
  }
  
};
