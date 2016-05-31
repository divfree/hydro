#include "linear_algebra.hpp"



bool check_validity_double(const double& a)
{
  return (a*0==0);
}

bool check_validity(const double& a)
{
  return check_validity_double(a);
}

bool hydro2D_uniform::check_validity(const double& a)
{
  return check_validity_double(a);
}


// retrun b : b=2^k-1, b>=a, k->min
unsigned int sup_power2(unsigned int a)
{
  unsigned int b=1;
  while(a>b-1)
  {
    b = b<<1;
  }
  return b;
}

double get_nan()
{
  double a=0.0;
  return 0.0/a;
}
