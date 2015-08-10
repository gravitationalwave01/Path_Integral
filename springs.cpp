//This code will solve the problem of computing the free energy change when a system goes from 

// |--o--o--o--| (lambda = 0) to |--o--o--o--------| (lambda = 1)
//where | is a wall, -- is a spring, and o is a classical block.
//thus we are studying how the free energy changes when we increaase the total spring length
//we will do this by computing the work done in this process and averging using Jarzynski

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


class Pos
{
public:
  double start, end;

  Pos(double s, double e)
  {
    start = s;
    end = e;
  }
  
  double get_pos(double lambda)
  {
    return start+(end-start)*lambda;
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
  double eq_length;

  RP(int d, int numBeads, std::vector<double>* positions, std::vector<double>* velocities, std::vector<double>* masses, double o, double eq)
  {
    P = numBeads;
    dim = d;
    omega = o;
    eq_length = eq;

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
int velocityVerlet(RP* myRP, double deltaT, double wall)
{
  std::vector<double> curForces;
  std::vector<double> oldForces;
  curForces.resize(myRP->P);
  oldForces.resize(myRP->P);
  double dim = myRP->dim;
  double tmp;
  double omega = myRP->omega;
  double eq_length = myRP->eq_length;

  //update forces:
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
    {
      if (curBead == 0)
	tmp = - omega*omega*myRP->beads[0].mass * ( (myRP->beads[0].position[j] - eq_length) - (myRP->beads[1].position[j] - myRP->beads[0].position[j] - eq_length) );
      else if (curBead == myRP->P - 1)
	tmp = - omega*omega*myRP->beads[myRP->P-1].mass * 
	        ( (myRP->beads[myRP->P-1].position[j] - myRP->beads[myRP->P-2].position[j] - eq_length) 
	        - (wall - myRP->beads[myRP->P-1].position[j] - eq_length ) );
      else 
	tmp = - omega*omega*myRP->beads[curBead].mass * 
	        ( (myRP->beads[curBead].position[j] - myRP->beads[curBead-1].position[j] - eq_length) 
	        - (myRP->beads[curBead+1].position[j] - myRP->beads[curBead].position[j] - eq_length) );
      
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
      if (curBead == 0)
	tmp = - omega*omega*myRP->beads[0].mass * ( (myRP->beads[0].position[j] - eq_length) - (myRP->beads[1].position[j] - myRP->beads[0].position[j] - eq_length) );
      else if (curBead == myRP->P - 1)
	tmp = - omega*omega*myRP->beads[myRP->P-1].mass * 
	        ( (myRP->beads[myRP->P-1].position[j] - myRP->beads[myRP->P-2].position[j] - eq_length) 
		- (wall - myRP->beads[myRP->P-1].position[j] - eq_length ) );
      else 
	tmp = - omega*omega*myRP->beads[curBead].mass * 
	        ( (myRP->beads[curBead].position[j] - myRP->beads[curBead-1].position[j] - eq_length) 
		- (myRP->beads[curBead+1].position[j] - myRP->beads[curBead].position[j] - eq_length) );
      
      curForces[curBead*dim+j] = tmp;
    }
  
  //update velocities
  for (int curBead = 0; curBead < myRP->P; curBead++)
    for (int j=0;j<dim;j++)
      myRP->beads[curBead].velocity[j] += deltaT * 0.5 / myRP->beads[curBead].mass * (oldForces[curBead*dim+j]+curForces[curBead*dim+j]);
  
};







int main()
{
  //initialize some constants
  int dim = 1;
  int P = 3; //number of blocks
  double timevec[6] = {1,2,3,5,10,20}; // total time to go from lambda = 0 to lambda = 1
  double time;
  double deltaT = 0.01;
  double collFreq = 0.01; //frequency of collision with bath
  double eq_length = 1; //equilibrium length of each spring
  double omega = 1; // spring stiffness
  double answer = 1;
  double P1[2] = {1,1};
  double P2[2] = {2,1};


  double ibeta = 0.1; // ibeta = 1/beta = kT
  //Initialize the random number generators
  std::default_random_engine generator(rdtsc());
  std::normal_distribution<double> pdist(0,sqrt(ibeta));
  std::uniform_real_distribution<double> unif(0.0,1.0);

  int numTrials = 6;
  std::vector<double> deltaF;

  //open the output files:
  std::ofstream fpos("positions.dat"); //for positions
  //  std::ofstream fvel("velocities.dat"); //for velocities
  //    std::ofstream fene("energy.dat"); //for energies (kinetic and potential)
  std::ofstream fpot("fpot.dat");   // for dU/dlambda
  std::ofstream dist("dist.dat"); //stores the distribution of velocities (should be an MB distribution)
  std::stringstream fname;
  /*
    if (trial % 2 == 0)
    fname << "fwork_t" << trial/2 << ".dat";
    else
    fname << "rwork_t" << trial/2 << ".dat";
  */
  fname << "work.dat";
  std::ofstream work(fname.str().c_str()); //stores the computed values for the work
  int obs_trial = 2; //the trial about which we record dyanmics, work, distributions, etc.

  for (int trial = 0; trial < numTrials; trial++)
  {
    time = timevec[trial];
    Pos pos((P+P1[0])*eq_length, (P+P2[0])*eq_length);
    std::cout << "starting trial " << trial+1 << " of " << numTrials << std::endl;



    //Initialize a vector that will store an ensemble of ring polymers
    int numPolymers = 10000;
    std::vector<RP> ensemble;
    double tmp = 0;


    //initialize the ensemble with RP velocities randomly from this distribution
    //this equilibrates the initial state
    for (int i = 0; i < numPolymers; i++)
    {
      std::vector<double> positions;
      std::vector<double> velocities;
      std::vector<double> masses;
    
      for (int curBead = 0; curBead < P; curBead++)
      {
	masses.push_back(1.0);
	for (int curDim = 0; curDim < dim; curDim++)
	{
	  velocities.push_back(pdist(generator));
	  positions.push_back((curBead+1)*eq_length); //equilibrium position
	}	
      }
    
      //create a polymer with these initial conditions and store in the ensemble 
      ensemble.push_back(RP(dim,P,&positions,&velocities,&masses,omega,eq_length));
    }
  
    //we have the ensemble initialized with the proper kinetic energy, but we still need to get the right configuration
    //for this, we propogate in time for a bit:
    for (double time_t = 0; time_t < 10; time_t += deltaT) //we don't care too much about energy conservation -- just get a distribution!
    {
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	velocityVerlet(&*myRP, deltaT, pos.get_pos(0));
    }
    
    for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
      if (trial == obs_trial)
	dist << myRP->beads[0].position[0] << std::endl;
      
    //now the kinetic energies are all messed up, so re-assign them:
    for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
      for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	for (int j=0;j<dim;j++)
	  it->velocity[j] = pdist(generator); //only while the mass = 1
    


    
    double pot;//pot stores values of dU/dlambda
    double dlambda = deltaT/time;
    double exp_work=0;
    double lambda = 0;
    
    std::vector<std::vector<std::pair<double, double> > > integrand; 
    for (int i=0;i<numPolymers;i++)
    {
      std::vector<std::pair<double,double> > blah;
      integrand.push_back(blah);
    } 
    std::vector<std::vector<std::pair<double, double> > >::iterator it_integrand;
    
    

    lambda = 0;
    for (double curTime = 0; curTime <= time; curTime+=deltaT)
    {
      
      if (trial == obs_trial)
	fpos << lambda << " " << ensemble[0].beads[0].position[0] << " " << ensemble[0].beads[1].position[0] 
	     << " " << ensemble[0].beads[2].position[0] << std::endl;
      
      
      //Andersen thermostat
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	for (std::vector<Bead>::iterator it = myRP->beads.begin(); it < myRP->beads.end(); it++)
	  if (unif(generator) < deltaT * collFreq) 
	    for (int j=0;j<dim;j++)
	      it->velocity[j] = pdist(generator); //only while the mass = 1
      
      it_integrand = integrand.begin();

      //propagate in time
      for(std::vector<RP>::iterator myRP = ensemble.begin(); myRP < ensemble.end(); myRP++)
	{

	  //int velocityVerlet(RP* myRP, double deltaT, double wall)
	  velocityVerlet(&*myRP,deltaT,pos.get_pos(lambda));
	
	  // sample dU/dlambda with some probability.
	  // to sample dU/dlambda we really sample x, plug this x into (U_{lambda+dlambda/2) - U_{lambda-dlambda}) / 2
	  pot = 0;
	  for (int j=0;j<dim;j++)
	    pot += pos.get_pos(lambda) - myRP->beads[P-1].position[j];
	  //	  pot += ( pow(myRP->beads[P-1].position[j]-pos.get_pos(lambda+dlambda),2) - pow(myRP->beads[P-1].position[j]-pos.get_pos(lambda-dlambda),2) );
	  pot *= myRP->beads[P-1].mass*omega*omega;
	
	//	std::cout << pot <<std::endl;
        it_integrand->push_back(std::pair<double,double>(curTime,pot)); 
        it_integrand++;

      } //end of ensemble loop
      
	lambda += dlambda;
    } // end of time loop
    

    tmp = 0;
    double tmp_avg=0;
    for (it_integrand = integrand.begin();it_integrand<integrand.end();it_integrand++)
    {
      tmp = 1/time * simpson38(&*it_integrand);
      if (trial == obs_trial)
	work << tmp << std::endl;
      exp_work += exp(-1/ibeta * tmp); //assumes that mass=1 therefore a^2 = kT
      tmp_avg += tmp;
    }

    std::cout << "simulation time is " << time << std::endl;
    std::cout << "W^a is " << tmp_avg / numPolymers << std::endl;
    std::cout << "W^z is " << -ibeta*log(exp_work/numPolymers) << std::endl;
      //tmp_out << tmp_avg / numPolymers << " " << -ibeta*log(exp_work/numPolymers) << std::endl;

    //close all output files:
    //  fpos.close();
    //  fvel.close();
    //  fene.close();
    //  fpot.close();

  }//end of trial loop
  

};
