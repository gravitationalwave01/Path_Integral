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

double getForce(double position, double mass)
{
  double omega = 1;
  return -mass*omega*omega*position;
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


/*
double combine(double pos, double lambda, HarmonicPotential* p1, HarmonicPotential* p2)
{
  return (1-lambda)*p1->getPotential(pos) + lambda*p2->getPotential(pos);
}
*/
/*
class TwoStatePotential
{
 public:
  double getPotential(double lambda, TwoStatePotential::*pot_A, TwoStatePotential::*pot_B, double x);
  void harmonic(double omega, double mass, double pos);
};

TwoStatePotential::harmonic(double omega, double mass, double pos)
{
  return 0.5*mass*omega*omega*pos*pos;
}

OBTwoStatePotential::getPotential(double lambda, TwoStatePotential::*pot_A, TwoStatePotential::*pot_B, double x)
{
  
}
*/

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

  /*
  void updateForces()
  {
    for (int i=0;i<dim;i++)
      force[i] = getForce(position[i],mass);
  }
  */

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
    //    std::cout << "numbeads is " << numBeads << std::endl;
    for (int i =0; i < numBeads; i++)
    {
      //      std::cout << "here" << std::endl;
      //      std::cout << *(positions+i*dim) << " " << dim << " " << *(velocities+i*dim) << " " << masses[i] << std::endl;
      //      Bead* b = new Bead(positions+i*dim, dim, velocities+i*dim, masses[i]);
      //      b->printMe();
      beads.push_back(Bead(positions->begin()+i*dim, dim, velocities->begin()+i*dim, masses->at(i)));
    }
    //    std::cout << "supposedly initialized all the beads" << std::endl;
    //    beads.at(0).printMe();
  }

  /*
  void updateForces()
  {
    for (std::vector<Bead>::iterator it = beads.begin(); it < beads.end(); it++)
      it->updateForces();
  }
  */
  
  void printMe()
  {
    std::cout << "DETAILS OF THE RING POLYMER::" << std::endl;
    std::cout << "POSITIONS:     VELOCITIES:" << std::endl;
    //    for (int i = 0; i < P; i++)
    //      std::cout << beads[i].position << "                       " << beads[i].velocity << std::endl;
  }
};


double integrate(std::vector<std::pair<double,double> >* func)
{
  std::cout << "starting integration" << std::endl;
  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly and constantly
  std::cout << delta << std::endl;
  std::vector<std::pair<double,double> > derivative;
  derivative.push_back(std::pair<double,double>(func->at(0).first, (func->at(1).second - func->at(0).second)/(delta)));

  //  for (std::vector<std::pair<double,double> >::iterator it = func->begin(); it+2 < func->end(); it++)
  //std::cout << it->first << " " << it->second << std::endl;

  double tmp = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin(); it+2 < func->end(); it++)
  {
    //    std::cout <<  ((it+2)->second - it->second)/(2*delta) << std::endl;
    derivative.push_back(std::pair<double,double>((it+1)->first, ((it+2)->second - it->second)/(2*delta)));
  }
  derivative.push_back(std::pair<double,double>((func->end()-1)->first, ((func->end()-1)->second - (func->end()-2)->second)/(delta)));

  //  std::cout << "second to last : " << (func->end()-1)->second << std::endl;
  //  std::cout << "last : " << func->end()->second << std::endl;

  //  std::cout << "This is what I got for the derivative : " << std::endl;
  //  for (std::vector<std::pair<double,double> >::iterator it = derivative.begin(); it < derivative.end(); it++)
  //    std::cout << it->first << " " << it->second << std::endl;
  //  std::cout << " end of derivative" << std::endl;

  std::cout << "computed the derivative. starting FE integration" << std::endl;
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = derivative.begin(); it+1 < derivative.end(); it++)
  {
    ans += delta * ( (it+1)->second + it->second) / 2;
    //    std::cout << "between " << it->first << " and " << (it+1)->first << " we get " << ans << std::endl;
  }
  
  return ans;
}


int main()
{
  //initialize the ring polymer 
  int dim = 1;
  int P = 3;
  //  double* positions = new double[dim*P];
  //  double* velocities = new double[dim*P];
  //  double* masses = new double[dim];
  //  positions[0] = 1;
  //  positions[1] = 2;
  //  positions[2] = 3;
  //  velocities[0] = 0;
  //  velocities[1] = 0;
  //  velocities[2] = 0;
  //  masses[0] = 1;

  //  RP myRP(dim,P,positions,velocities,masses);
  double time = 300;
  double deltaT = 0.02;

  //open the output files:
  std::ofstream fpos("positions.dat");
  std::ofstream fvel("velocities.dat");
  std::ofstream fene("energy.dat");
  

  //create the two potentials by first initializing their classes:
  HarmonicPotential hp1(1,2);
  HarmonicPotential hp2(1,3);

  /*
  //create function pointers to the two potential functions:
  double (*pot_A)(double);
  pot_A = &hp1.getPotential;
  double (*pot_B)(double);
  pot_B = &hp2.getPotential;

  //create function pointers to the two force functions:
  double (*force_A)(double);
  force_A = &hp1.getForce;
  double (*force_B)(double);
  force_B = &hp2.getForce;
  */  

  /*
  std::cout << "TESTING FUNCTION OF FUNCTION FUNCTIONALITY" << std::endl;
  for (double lambda = 0; lambda <= 1; lambda+=0.1)
  {
    std::cout << combine(1,lambda, &hp1, &hp2) << std::endl;
  }
  std::cout << "DONE TESTING FUNCTION OF FUNCTION FUNCTIONALITY" << std::endl;
  */

  //For each instance of the parameter, create a boltzmann distribution of polymers:
  int numPolymers = 1000;
  std::vector<RP> ensemble;
  //create a boltzmann distribution random number generator:
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0,0.25);
  
  //initialize many polymers with positions and velocities randomly from this distribution
  //make sure that the center of mass (COM) momentum of each polymer is 0
  double com = 0;
  std::cout << " **** STARTING RP INIT  ****" << std::endl;
  for (int i = 0; i < numPolymers; i++)
  {
    com = 0;
    std::vector<double> positions;
    std::vector<double> velocities;
    std::vector<double> masses;
    
    for (int curBead = 0; curBead < P; curBead++)
    {
      masses.push_back(1);
      for (int curDim = 0; curDim < dim; curDim++)
      {
	positions.push_back(distribution(generator));
	double cur_v = distribution(generator);
	velocities.push_back(cur_v);
	com += cur_v; //WARNING: ASSUMES THAT THE MASS IS 1!!!!!!!!!!!
      }	
    }
    
    //subtract the COM momentum:
    for (int curBead = 0; curBead < P; curBead++)
      for (int curDim = 0; curDim < dim; curDim++)
	velocities[curBead*dim+curDim] -= com/(dim*P); //ONLY WORKS WHEN MASS IS 1!!!!!!!!!!!!
    
    //create a polymer with these initial conditions and store in the ensemble 
    ensemble.push_back(RP(dim,P,&positions,&velocities,&masses));
  }
  
  
  //output the initial distribution of velocities and positions:
  std::ofstream dist("dist.dat");
  for (int curRP = 0; curRP < numPolymers; curRP++)
  {
    for (int curBead = 0; curBead < P; curBead++)
    {
      for (int curDim = 0; curDim < dim; curDim++)
      {
	dist << ensemble[curRP].beads[curBead].position[curDim] << " " << ensemble[curRP].beads[curBead].velocity[curDim] << std::endl;
      }	
    }
  }
  dist.close();
  
  //  srand (time(NULL));
  
  //we have all of the ring polymers initialized -- now propogate them in time so that they equilibrate.
  std::vector<double> curForces; 
  std::vector<double> oldForces; 
  double kinetic, potential;
  bool first = true;
  double tmp;
  std::ofstream fpot("fpot.dat");
  
  double energy_test = 0;
  bool testing = true;
  std::vector<std::pair<double, double> > integrand; //will store the values (lambda, potential energy)
  

  //as we vary lambda, the potential changes and we need to wait for the system to equlibrate. 
  //after equlibrating we can sample the potential energy
  int numSamples = 0;
  double pot = 0;

  // Smoothly vary the potential-flipping parameter:
  for (double lambda = 0; lambda <= 1; lambda += 0.02)
  {
    pot = 0;
    numSamples = 0;
    std::cout << " ************** CURRENT LAMBDA ************* " << lambda << std::endl;    
    for (double curTime = 0; curTime < time; curTime+=deltaT)
    {
      kinetic=0;
      potential=0;
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
	for (int curBead = 0; curBead < myRP->P; curBead++)
	  for (int j=0;j<dim;j++)
	  {
	    kinetic += 0.5 * myRP->beads[curBead].mass * myRP->beads[curBead].velocity[j] * myRP->beads[curBead].velocity[j];
	    potential += ( (1-lambda)*hp1.getPotential(myRP->beads[curBead].position[j]) + lambda*hp2.getPotential(myRP->beads[curBead].position[j]) );
	  }



	/*
	if (testing)
	{
	  energy_test = 0;
	  //	  if (lambda <0.05)
	  {
	    energy_test += 0.5 * ensemble[0].beads[0].mass * ensemble[0].beads[0].velocity[0] * ensemble[0].beads[0].velocity[0];
	    energy_test += ( (1-lambda)*hp1.getPotential(ensemble[0].beads[0].position[0]) + lambda*hp1.getPotential(ensemble[0].beads[0].position[0]) );
	    std::cout << curTime + lambda*1000 << " " << energy_test << " " << ensemble[0].beads[0].position[0] << " " <<ensemble[0].beads[0].velocity[0] << std::endl;
	  }
	  testing = false;
	}
	*/
      } //end of VV loop
      
      fpos << curTime+time*lambda*50 << " ";
      fvel << curTime+time*lambda*50 << " ";

      //record the position and velocity of just one of the polymers:
      for (std::vector<Bead>::iterator it = ensemble[1].beads.begin(); it < ensemble[1].beads.end(); it++)
      {
	for (int j=0;j<dim;j++)
	{
	  fpos << it->position[j] << " "; 
	  fvel << it->velocity[j] << " ";
	}
      }
      fpos << std::endl;
      fvel << std::endl;
      fene << curTime+time*lambda*50 << " " << kinetic << " " << potential << std::endl; 

      // sample the potential energy with some probability:
      if (distribution(generator) >= 0.1 && curTime > 1) // 
      {
	pot += potential;
	numSamples++;
      }

    } //end of time loop
    integrand.push_back(std::pair<double,double>(lambda,pot));

    fpot << lambda << " " << pot/numSamples << " " << numSamples << std::endl;
  } // end of lambda loop
  fpos.close();
  fvel.close();
  fene.close();
  fpot.close();
  std::cout << "FREE ENERGY DIFFERENCE IS : " << std::endl;

  for (std::vector<std::pair<double,double> >::iterator it = integrand.begin(); it < integrand.end(); it++)
    std::cout << it->first << " " << it->second << std::endl;


  std::vector<std::pair<double,double> > test;
  /*
  test.push_back(std::pair<double,double>(0.0,0));
  test.push_back(std::pair<double,double>(0.5,.25));
  test.push_back(std::pair<double,double>(1.0,1));
  test.push_back(std::pair<double,double>(1.5,2.25));
  test.push_back(std::pair<double,double>(2.0,4));
  test.push_back(std::pair<double,double>(2.5,6.25));
  test.push_back(std::pair<double,double>(3.0,9));
  test.push_back(std::pair<double,double>(3.5,12.25));
  test.push_back(std::pair<double,double>(4.0,16));
  test.push_back(std::pair<double,double>(4.5,20.25));
  test.push_back(std::pair<double,double>(5.0,25));
  */  
  std::cout << integrate(&integrand) << std::endl;
  
};
