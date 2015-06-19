#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <math.h>

int main()
{
  /*
  std::default_random_engine generator;
  double a = 1; // a^2 = kT/m
  std::normal_distribution<double> distribution(0,a);

  std::ofstream output("test.dat");
  int numSample = 1000000;
  for (int i=0;i<numSample;i++)
    output << distribution(generator) << std::endl;
  */

  std::ofstream prob("test.dat");

  const double nrolls=1000000;  // number of experiments
  
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,.5);

  double p[10000]={};
  
  int counter = 0;
  for (int i=0; i<nrolls; ++i) {
    double number = distribution(generator);
    //    int num = (int)floor(number*100)+500;
    //    if (num < 0 || num > 1000) std::cout << "out of bounds" << std::endl;
    if ((number>=-5.0)&&(number<5.0)) 
    {
      p[(int)floor(number*1000)+5000]++; 
      counter++;
    }
    else std::cout << "out of bounds" << std::endl;
    
  }
  std::cout << "counter - numrolls = " << counter-nrolls << std::endl;
  //  std::cout << "normal_distribution (5.0,2.0):" << std::endl;
  
  for (int i=0; i<10000; i++) {
    //    std::cout << i << "-" << (i+1) << ": ";
    prob << ((double)i-5000.0)/1000 << " " <<  p[i]/nrolls << std::endl;
  }
  
  return 0;
}
