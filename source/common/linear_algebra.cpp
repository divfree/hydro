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


matrix<double> HilbertMatrix(int n)
{
  matrix<double> res(n,n);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
      res[i][j]=1.0/(i+j+1);
  return res;
}

column<double> RandomColumn(int n)
{
  //srand((unsigned int)time(NULL));
  column<double> res(n);
  for(int i=0; i<res.N; i++)
  {
    res[i]=(rand()%2000)/1000.0-1.0;
  }
  return res;
}

matrix<double> RandomMatrix(int m, int n)
{
  //srand((unsigned int)std::time(NULL));
  matrix<double> res(m,n);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
    {
      res[i][j]=(rand()%2000)/1000.0-1.0;
    }
  return res;
}


matrix<double> RandomSimilarMatrix(matrix<double> B)
{
  int N=B.n;
  matrix<double> C=RandomMatrix(N,N);
  matrix<double> Q,R;
  Householder(C, Q, R);
  // Q - some orthogonal matrix
  return tp(Q)*B*Q; // some matrix with eigenvalues from D
}

matrix<double> Sample1(int n)
{
  matrix<double> res(n,n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    {
      res[i][j]=max(i,j)+1;
    }
  return res;
}

matrix<double> Sample8(int n)
{
  matrix<double> res(n,n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    {
      res[i][j]=(n-max(i,j)+1);
    }
  return res;
}

matrix<double> Sample5(int n)
{
  matrix<double> res(n,n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    {
      res[i][j]=(i==0)?1:1.0/(i+j+1);
    }
  return res;
}

double energetic_norm(const matrix<double>& A, const column<double>& x)
{
  return sqrt(tp(A*x)*x);
}

matrix<double> zero(int m, int n)
{
  matrix<double> res(m,n);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
    {
      res[i][j]=0;
    }
  return res;
}

column<double> zero(int n)
{
  column<double> res(n);
  for(int i=0; i<res.N; i++)
  {
    res[i]=0;
  }
  return res;
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

/**********************************************************\
  Divide segment [a,b] into equal parts
  N - number of subintervals
\**********************************************************/
column<double> getGridRegular(double a, double b, int N)
{
  column<double> X=column<double>(N+1);
  for(int i=0; i<=N; i++) X[i]=a+(b-a)*((double)i/N);
  return X;
}

double get_nan()
{
  double a=0.0;
  return 0.0/a;
}
