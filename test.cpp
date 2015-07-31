#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <math.h>


double simpson(std::vector<std::pair<double,double> >* func) // implementing simpson's rule to integrate func
{
  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly 
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin() + 2; it < func->end(); it+=2)
  {
    ans += delta * ( (it-2)->second + 4 * (it-1)->second + it->second ) / 3; //simpson's 3/8 rule -- integrates three intervals at a time
    std::cout << ans << std::endl;
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
    std::cout << ans << std::endl;
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




//this function integrates *func using the traezoidal method.
double integrate(std::vector<std::pair<double,double> >* func)
{
  double delta = func->at(1).first - func->at(0).first; //assuming that the domain varies linearly and constantly
  double ans = 0;
  for (std::vector<std::pair<double,double> >::iterator it = func->begin(); it+1 < func->end(); it++)
    ans += delta * ( (it+1)->second + it->second) / 2;
  
  return ans;
}


unsigned long long rdtsc()
{
  unsigned int lo,hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return ((unsigned long long)hi << 32) | lo;
}


int main()
{
  std::vector<std::pair<double,double> > test;

  test.push_back(std::pair<double,double>(0.0,5));
  test.push_back(std::pair<double,double>(0.5,5.25));
  test.push_back(std::pair<double,double>(1.0,6));
  test.push_back(std::pair<double,double>(1.5,7.25));
  test.push_back(std::pair<double,double>(2.0,9));
  test.push_back(std::pair<double,double>(2.5,11.25));
  test.push_back(std::pair<double,double>(3.0,14));
  test.push_back(std::pair<double,double>(3.5,17.25));
  test.push_back(std::pair<double,double>(4.0,21));
  test.push_back(std::pair<double,double>(4.5,25.25));
  test.push_back(std::pair<double,double>(5.0,30));
  test.push_back(std::pair<double,double>(5.5,35.25));

  std::cout << "FREE ENERGY DIFFERENCE IS : " << std::endl;
  std::cout << "simpson's 3/8: " << simpson38(&test) << std::endl;
  std::cout << "simpson's: " << simpson(&test) << std::endl;
  std::cout << "trapezoid: " << integrate(&test) << std::endl;


  std::ofstream output("test.dat");
  double ibeta = 1.5;
  std::default_random_engine generator(rdtsc());
  std::normal_distribution<double> distribution(0,sqrt(ibeta));
  for (int i=0;i<1000000;i++)
    output << distribution(generator) << std::endl;
  
  return 0;
}
