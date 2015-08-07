// computing free energies differences of HO using MC


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
  int MC_steps_vec[3] = {10,100,500};
  int MC_steps;

  //  double a = 1.2247; // a^2 = kT/m
  double ibeta = 1.5; //inverse beta, ibeta = kT

  //create the two potentials
  double omega1 = 1;
  double omega2 = 2;
  HarmonicPotential hp1(1.0,omega1);
  HarmonicPotential hp2(1.0,omega2);
  double answer = ibeta * log(omega2/omega1);
  std::ofstream tmp_out("tmp.dat");

  int numTrials = 3;
  std::vector<double> deltaF;
  for (int trial = 0; trial < numTrials; trial++)
  {
    std::cout << "starting trial " << trial << " of " << numTrials << std::endl;
    MC_steps = MC_steps_vec[trial];
    std::cout << "MC_steps is " << MC_steps << std::endl;
    tmp_out << MC_steps << " ";
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
  std::uniform_real_distribution<double> unif(-1.0,1.0);
  std::uniform_real_distribution<double> unif01(0,1.0);


  //Initialize a vector that will store an ensemble of ring polymers
  int numPolymers = 100000;
  std::vector<RP> ensemble;


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
  
  
  double pot;//pot stores dU/dlambda


  //as we vary lambda, the potential changes 
  double dlambda = 1.0/MC_steps;
  std::vector<double> W;
  for (int i=0;i<numPolymers;i++)
    W.push_back(0);
  std::vector<double>::iterator W_it = W.begin();

  std::vector<double> expW;
  for (int i=0;i<numPolymers;i++)
    expW.push_back(exp(0));
  std::vector<double>::iterator expW_it = expW.begin();

  double step = 0.1;
  double new_pos = 0;
  double deltaV = 0; //change in energy imposed by metropolis displacement

  for (double lambda = dlambda; lambda <= 1.0001; lambda+=dlambda)
  {    
    //    std::cout << "lambda is " << lambda << std::endl;
    
    //implement metropolis method:
    //NOTE: THIS CODE ONLY WORKS FOR 1D HO !!!!!!!! <-----------------
    // (doesn't generalize to ring polymer of any sort!!)
    for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
    {
      new_pos = unif(generator) * step + myRP->beads[0].position[0];
      deltaV = lambda * (hp2.getPotential(new_pos)-hp2.getPotential(myRP->beads[0].position[0])) + (1-lambda)*(hp1.getPotential(new_pos)-hp1.getPotential(myRP->beads[0].position[0]));
      if (deltaV < 0)
	myRP->beads[0].position[0] = new_pos;
      else if (exp(-1/ibeta * deltaV) < unif01(generator))
	myRP->beads[0].position[0] = new_pos;
    }
    
	
    //compute the work done by making the shift above.
    W_it = W.begin();
    for (std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
    {
      *W_it += dlambda * (hp2.getPotential(myRP->beads[0].position[0]) - hp1.getPotential(myRP->beads[0].position[0]) ) ;
      //      *W_it += hp2.getPotential(myRP->beads[0].position[0]) - hp1.getPotential(myRP->beads[0].position[0]) ;
      W_it++;
    }

    
  } // end of time loop

  //close all output files:
  //  fpos.close();
  //  fvel.close();
  fene.close();
  //  fpot.close();

  double tmp = 0;
  double exp_tmp = 0;
  for (int i=0;i<numPolymers;i++)
  {
    tmp += W[i];
    exp_tmp += -ibeta * log(1/ibeta * W[i]);
  }
  

  std::cout << "W^a is " << tmp / numPolymers << std::endl;
  std::cout << "W^x is " << exp_tmp / numPolymers << std::endl;


  /*
  double tmp_avg=0;
  for (it_integrand = integrand.begin();it_integrand<integrand.end();it_integrand++)
  {
    tmp = 1/time * simpson38(&*it_integrand);
    work << tmp << std::endl;
    exp_work += exp(-1/ibeta * tmp); //assumes that mass=1 therefore a^2 = kT
    tmp_avg += tmp;
  }
  */

  //  std::cout << "W^a is " << tmp_avg / numPolymers << std::endl;
  //  std::cout << "W^z is " << -ibeta*log(exp_work/numPolymers) << std::endl;
  //  tmp_out << tmp_avg / numPolymers << " " << -ibeta*log(exp_work/numPolymers) << std::endl;

  //  deltaF.push_back(-ibeta*log(exp_work/numPolymers));
  //  std::cout << "Ratio of collisions with the bath : "<< (1 / collFreq) / (time/deltaT) << std::endl;
  }//end of trial loop

  /*
  double finalanswer = 0;
  for (int i=0;i<numTrials;i++)
    finalanswer += deltaF[i];

  finalanswer/=numTrials;

  double stdDev = 0;
  for (int i=0;i<numTrials;i++)
    stdDev += (deltaF[i]-finalanswer)*(deltaF[i]-finalanswer);

  stdDev /=numTrials;
  stdDev = sqrt(stdDev);
  */
  
  //compute and print the free energy difference
  /*
  std::cout << "COMPUTED FREE ENERGY DIFFERENCE via " << numTrials << " trials IS : " << finalanswer << std::endl;
  std::cout << "EXPECTED " << answer << std::endl;
  std::cout << "ERROR " << ((finalanswer-answer)/answer)*100 << "%" << std::endl;
  std::cout << "STD DEV of delta F values is " << stdDev << std::endl;
  */
};
