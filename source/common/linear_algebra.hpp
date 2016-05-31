#pragma once

//#define _DEBUG // Grid range check
//#define _DEBUG2 // fluxes calculation info
//#define _DEBUG3 // TEquation_unified.sort_refs()

#include "../common/std_ref.hpp"
#include <cmath>
#include <cassert>
using std::abs;

double get_nan();

bool check_validity_double(const double& a);

// return b : b=2^k-1, b>=a, k->min
unsigned int sup_power2(unsigned int a);
bool check_validity(const double& a);

template<class T>
class column
{
public:
  T* D;
  int N;
  int vN; // virtual N (not real size)
  column()
  {
    N=0;
    vN=0;
    D=0;
  }
  column(int sN)
  {
    N=sN;
    D=new T[N];
  }
  column(const column<T>& source)
  {
    N=source.N;
    D=new T[N];
    for(int i=0; i<N; i++) D[i]=source.D[i];
  }
  void change_N(int newN)
  {
    if(newN==N) return;
    T* newD=new T[newN];
    for(int i=0; i<min(N,newN); i++) newD[i]=D[i];
    delete [] D;
    D=newD;
    N=newN;
  }
  ~column()
  {
    delete[] D;
  }
  T& operator[](int i)
  {
#ifdef __RANGE_CHECK
    assert(i >=0 && i < N);
#endif
    return D[i];
  }
  const T& operator[](int i) const
  {
#ifdef __RANGE_CHECK
    assert(i >=0 && i < N);
#endif
    return D[i];
  }
  const T& operator()(int i) const
  {
    return D[i];
  }
  column<T>& operator=(const column<T>& S)
  {
    if(S.N!=N)
    {
      delete[] D;
      N=S.N;
      D=new T[N];
    }
    for(int i=0; i<N; i++)
    {
      D[i]=S.D[i];
    }
    return *this;
  }
  bool operator==(const column<T>& S)
  {
    if(S.N!=N) return false;
    bool res=true;
    int j=0;
    while(res && j<N)
    {
      res=res && (D[j]==S.D[j]);
      j++;
    }
    return res;
  }
  bool operator!=(const column<T>& S)
  {
    return !operator==(S);
  }
  column<T> operator+(const column<T>& S)
  {
    column<T> res(max(N, S.N));
    for(int i=0; i<min(N,S.N); i++)
    {
      res.D[i]=D[i]+S.D[i];
    }
    return res;
  }
  column<T> operator-(const column<T>& S)
  {
    column<T> res(max(N, S.N));
    for(int i=0; i<min(N,S.N); i++)

    {
      res.D[i]=D[i]-S.D[i];
    }
    return res;
  }
  column<T> operator*(double a) const
  {
    column<T> res(N);
    for(int i=0; i<N; i++)
    {
      res.D[i]=D[i]*a;
    }
    return res;
  }
  column<T> operator/(double a)
  {
    column<T> res(N);
    for(int i=0; i<N; i++)
    {
      res.D[i]=D[i]/a;
    }
    return res;
  }
  T norm()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=D[i]*D[i];
    return sqrt(res);
  }
  T sqrnorm()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=D[i]*D[i];
    return res;
  }
  T norm1()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=abs(D[i]);
    return res;
  }
  T Cheb()
  {
    T res=0;
	  for(int i=0; i<N; i++) res=max(res,abs(D[i]));
    return res;
  }
  void insert_to_end(T a)
  {
    change_N(N+1);
    D[N-1]=a;
  }
};

template<class T>
class row
{
  T* D;
public:
  int N;
  row()
  {
    N=0;
    D=0;
  }
  row(int sN)
  {
    N=sN;
    D=new T[N];
  }
  row(row<T>& source)
  {
    N=source.N;
    D=new T[N];
    for(int i=0; i<N; i++) D[i]=source.D[i];
  }
  ~row()
  {
    delete[] D;
  }
  T& operator[](int i)
  {
    return D[i];
  }
  const T& operator()(int i) const
  {
    return D[i];
  }
  row<T>& operator=(const row<T>& S)
  {
    if(S.N!=N)
    {
      delete[] D;
      N=S.N;
      D=new T[N];
    }
    for(int i=0; i<N; i++)
    {
      D[i]=S.D[i];
    }
    return *this;
  }
  row<T> operator+(row<T>& S)
  {
    row<T> res(max(N, S.N));
    for(int i=0; i<min(N,S.N); i++)
    {
      res.D[i]=D[i]+S.D[i];
    }
    return res;
  }
  row<T> operator-(row<T>& S)
  {
    row<T> res(max(N, S.N));
    for(int i=0; i<min(N,S.N); i++)
    {
      res.D[i]=D[i]-S.D[i];
    }
    return res;
  }
  row<T> operator*(T a) const
  {
    row<T> res(N);
    for(int i=0; i<N; i++)
    {
      res.D[i]=a*D[i];
    }
    return res;
  }
  row<T> operator/(T a)
  {
    row<T> res(N);
    for(int i=0; i<N; i++)
    {
      res.D[i]=D[i]/a;
    }
    return res;
  }
  T norm()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=D[i]*D[i];
    return sqrt(res);
  }
  T Cheb()
  {
    T res=0;
	for(int i=0; i<N; i++) res=max(res,abs(D[i]));
    return res;
  }
};

template<class T>
row<T> tp(const column<T>& S)
{
  row<T> res(S.N);
  for(int i=0; i<S.N; i++) res[i]=S(i);
  return res;
}
template<class T>
column<T> tp(const row<T>& S)
{
  column<T> res(S.N);
  for(int i=0; i<S.N; i++) res[i]=S(i);
  return res;
}
template<class T>
row<T> operator*(T a, const row<T>& B)
{
  return B*a;
}
template<class T>
column<T> operator*(double a, const column<T>& B)
{
  return B*a;
}

template<class T>
class matrix
{
  T** M; // col<row<T>>
public:
  int m,n;
  matrix()
  {
    m=0; n=0;
    M=0;
  }
  matrix(int sm, int sn)
  {
    m=sm; n=sn;
    M=new T*[m];
    for(int i=0; i<m; i++) M[i]=new T[n];
  }
  matrix(const matrix<T>& S)
  {
    m=S.m; n=S.n;
    M=new T*[m];
    for(int i=0; i<m; i++) M[i]=new T[n];
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++)
        M[i][j]=S.M[i][j];
  }
  ~matrix()
  {
    for(int i=0; i<m; i++) delete [] M[i];
    delete [] M;
  }
  T*& operator[](int i)
  {
    return M[i];
  }
  const T& operator()(int i, int j) const
  {
    return M[i][j];
  }
  void swap_row(int i1, int i2)
  {
    for(int j=0; j<n; j++) swap(M[i1][j], M[i2][j]);
  }
  void swap_col(int j1, int j2)
  {
    for(int i=0; i<m; i++) swap(M[i][j1], M[i][j2]);
  }
  void plus_row(int i1, T ci1, int i2)
  {
    for(int j=0; j<n; j++) M[i2][j]+=M[i1][j]*ci1;
  }
  void plus_col(int j1, T cj1, int j2)
  {
    for(int i=0; i<m; i++) M[i][j2]+=M[i][j1]*cj1;
  }
  row<T> get_row(int i)
  {
    row<T> res(n);
    for(int j=0; j<n; j++)
    {
      res[j]=M[i][j];
    }
    return res;
  }
  column<T> get_col(int j)
  {
    column<T> res(m);
    for(int i=0; i<m; i++)
    {
      res[i]=M[i][j];
    }
    return res;
  }
  column<T> get_col(int j, int i1, int i2)
  {
    if(i1>i2) swap(i1,i2);
    column<T> res(i2-i1);
    for(int i=0; i<m; i++)
    {
      res[i]=M[i+i1][j];
    }
    return res;
  }
  column<T> get_diag()
  {
    column<T> res(n);
    for(int i=0; i<n; i++) res[i]=M[i][i];
    return res;
  }
  matrix<T>& operator=(const matrix<T>& S)
  {
    if(S.m!=m || S.n!=n)
    {
      for(int i=0; i<m; i++) delete [] M[i];
      delete [] M;
      m=S.m; n=S.n;
      M=new T*[m];
      for(int i=0; i<m; i++) M[i]=new T[n];
    }
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++)
      {
        M[i][j]=S(i,j);
      }
    return *this;
  }
  matrix<T> operator+(const matrix<T>& S)
  {
    matrix<T> res(max(m,S.m),max(n,S.n));
    for(int i=0; i<min(m,S.m); i++)
      for(int j=0; j<min(n,S.n); j++)
      {
        res[i][j]=M[i][j]+S(i,j);
      }
    return res;
  }
  matrix<T> operator-(const matrix<T>& S)
  {
    matrix<T> res(max(m,S.m),max(n,S.n));
    for(int i=0; i<min(m,S.m); i++)
      for(int j=0; j<min(n,S.n); j++)
      {
        res[i][j]=M[i][j]-S(i,j);
      }
    return res;
  }
  matrix<T> operator*(T a) const
  {
    matrix<T> res(m,n);
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++)
      {
        res[i][j]=M[i][j]*a;
      }
    return res;
  }
  matrix<T> operator/(T a)
  {
    matrix<T> res(m,n);
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++)
      {
        res[i][j]=M[i][j]/a;
      }
    return res;
  }
  matrix<T> operator*(const matrix<T>& S) const
  {
    if(n!=S.m) throw string("incompatible matrix sizes");
    matrix<T> res(m, S.n);
    for(int i=0; i<res.m; i++)
      for(int j=0; j<res.n; j++)
      {
        res[i][j]=0;
        for(int k=0; k<S.m; k++)
        {
          res[i][j]+=M[i][k]*S.M[k][j];
        }
      }
    return res;
  }
  column<T> operator*(const column<T>& x) const
  {
    if(n!=x.N) throw string("incompatible matrix sizes");
    column<T> y(m);
    for(int i=0; i<m; i++)
    {
      y[i]=0;
      for(int j=0; j<m; j++)
      {
        y[i]+=M[i][j]*x(j);
      }
    }
    return y;
  }
  matrix<T> operator*(const row<T>& x) const
  {
    if(n!=1) throw string("incompatible matrix sizes");
    matrix<T> res(m,x.N);
    for(int i=0; i<m; i++)
      for(int j=0; j<x.N; j++)
      {
        res[i][j]=M[i][0]*x[j];
      }
    return res;
  }
  matrix<T> null_small(T eps)
  {
    matrix<T> res(m,n);
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++)
      {
        if(abs(M[i][j])<eps) res[i][j]=0; else res[i][j]=M[i][j];
      }
    return res;
  }
  T mul_diag()
  {
    if(n!=m) throw string("non-square matrix");
    T t=1;
    for(int i=0; i<m; i++) t*=M[i][i];
    return t;
  }
  T trace()
  {
    if(n!=m) throw string("non-square matrix");
    T t=0;
    for(int i=0; i<m; i++) t+=M[i][i];
    return t;
  }
	T Cheb()
	{
		T smax=0;
		for(int i=0; i<m; i++)
		{
			T sum=0;
			for(int j=0; j<n; j++)
			{
				sum+=abs(M[i][j]);
			}
			smax=max(smax,sum);
		}
		return smax;
	}
};

namespace hydro2D_uniform
{
bool check_validity(const double& a);

template<class T>
class grid
{
  T* M; // [x][y]
public:
  int Nx,Ny;
  grid()
  {
    Nx=0; Ny=0;
    M=0;
  }
  grid(int sNx, int sNy)
  {
    Nx=sNx; Ny=sNy;
    M=new T[Nx*Ny];
  }
  grid(const grid<T>& S)
  {
    Nx=S.Nx; Ny=S.Ny;
    M=new T[Nx*Ny];
    for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
    {
      M[Ny*i+j]=S(i,j);
    }
  }
  ~grid()
  {
    delete [] M;
  }
  T& operator()(int i, int j)
  {
    #ifdef _DEBUG
      if(i<0 || i>=Nx || j<0 || j>=Ny) throw string("Access violation: (i,j)="+point_to_string(i,j)+" (Nx,Ny)="+point_to_string(Nx,Ny));
    #endif
    return M[Ny*i+j];
  }
  const T& operator()(int i, int j) const
  {
    return M[Ny*i+j];
  }
  grid<T>& operator=(const grid<T>& S)
  {
    if(S.Nx!=Nx || S.Ny!=Ny)
    {
      Nx=S.Nx; Ny=S.Ny;
      delete [] M;
      M=new T[Nx*Ny];
    }
    for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
    {
      (*this)(i,j)=S(i,j);
    }
    return *this;
  }
  void change_N(int sNx, int sNy)
  {
    if(sNx!=Nx || sNy!=Ny)
    {
      operator=(grid<T>(sNx,sNy));
    }
  }
  grid<T> operator+(const grid<T>& S) const
  {
    if(S.Nx!=Nx || S.Ny!=Ny) throw string("Incompatible sizes of grids for operator+");
    grid<T> res(Nx,Ny);
    for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
    {
      res(i,j)=(*this)(i,j)+S(i,j);
    }
    return res;
  }
  grid<T> operator-(const grid<T>& S) const
  {
    if(S.Nx!=Nx || S.Ny!=Ny) throw string("Incompatible sizes of grids for operator-");
    grid<T> res(Nx,Ny);
    for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
    {
      res(i,j)=(*this)(i,j)-S(i,j);
    }
    return res;
  }
  grid<T> operator*(double a) const
  {
    grid<T> res(Nx,Ny);
    for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
    {
      res(i,j)=(*this)(i,j) * a;
    }
    return res;
  }
  grid<T> operator/(double a) const
  {
    grid<T> res(Nx,Ny);
    for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
    {
      res(i,j)=(*this)(i,j) / a;
    }
    return res;
  }
};

template<class T>
ostream& operator<<(ostream& out, grid<T> M)
{
  for(int j=M.Ny-1; j>=0; j--)
  {
    for(int i=0; i<M.Nx; i++)
    {
      out<<M(i,j)<<" ";
    }
    out<<endl;
  }
  return out;
}

template<class T>
bool check_validity(const grid<T>& G)
{
  bool res=true;
  for(int i=0; i<G.Nx; i++)
  for(int j=0; j<G.Ny; j++)
  {
    res = res && check_validity(G(i,j));
  }
  return res;
}

template<class T>
bool check_validity(const grid<T>& G, int i1, int j1, int i2, int j2)
{
  bool res=true;
  for(int i=i1; i<=i2; i++)
  for(int j=j1; j<=j2; j++)
  {
    res = res && check_validity(G(i,j));
  }
  return res;
}

// namespace end
}

template<class T>
matrix<T> tp(const matrix<T>& S)
{
  matrix<T> res(S.n, S.m);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
    {
      res[i][j]=S(j,i);
    }
  return res;
}
template<class T>
row<T> operator*(const row<T>& x, const matrix<T>& A)
{
  if(x.N!=A.m) throw string("incompatible matrix sizes");
  row<T> y(A.n);
  for(int j=0; j<y.N; j++)
  {
    y[j]=0;
    for(int i=0; i<A.m; i++)
    {
      y[j]+=x(i)*A(i,j);
    }
  }
  return y;
}
template<class T>
matrix<T> operator*(const column<T>& x, const matrix<T>& A)
{
  if(A.m!=1) throw string("incompatible matrix sizes");
  matrix<T> res(x.N, A.n);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
    {
      res[i][j]=x(i)*A(0,j);
    }
  return res;
}
template<class T>
T operator*(const row<T>& x, const column<T>& y)
{
  if(x.N!=y.N) throw string("incompatible matrix sizes");
  T res=0;
  for(int i=0; i<x.N; i++)
  {
    res+=x(i)*y(i);
  }
  return res;
}

template<class T>
matrix<T> operator*(column<T>& x, row<T>& y)
{
  matrix<T> res(x.N, y.N);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
    {
      res[i][j]=x[i]*y[j];
    }
  return res;
}

template<class T>
matrix<T> operator*(double t, const matrix<T>& A)
{
  matrix<T> res(A.m, A.n);
  for(int i=0; i<res.m; i++)
    for(int j=0; j<res.n; j++)
    {
      res[i][j]=t*A(i,j);
    }
  return res;
}

/*template<class T>
column<T> operator*(double t, column<T> V)
{
  column<T> res(V.N);
  for(int i=0; i<res.N; i++)
  {
    res[i]=t*V[i];
  }
  return res;
}*/

template<class T>
ostream& operator<<(ostream& out, matrix<T> M)
{
  for(int i=0; i<M.m; i++)
  {
    for(int j=0; j<M.n; j++)
    {
      out<<M[i][j]<<" ";
    }
    out<<endl;
  }
  return out;
}

template<class T>
istream& operator>>(istream& in, column<T>& C)
{
  char p;
  in>>p;
  if(p!='(') throw string("Unexpected symbol '")+p+string("' instead of '('");

  C=column<T>();
  int n=0;
  if(in.peek()!=')')
  {
    while(p!=')')
    {
      T a;
      in>>a;
      if(in.fail()) throw string("Reading of element failed");
      if(n>=C.N) C.change_N(sup_power2(n+1));
      C[n]=a;
      in>>p;
      n++;
    }
    C.change_N(n);
  }
  else in>>p;
  return in;
}

template<class T>
ostream& operator<<(ostream& out, column<T> C)
{
  out<<"(";
  for(int i=0; i<C.N; i++)
  {
      out<<(i==0?"":",")<<C[i];
  }
  out<<")";
  return out;
}


template<class T>
ostream& operator<<(ostream& out, row<T> R)
{
  for(int i=0; i<R.N; i++)
  {
      out<<R[i]<<" ";
  }
  return out;
}

template<class T>
matrix<T> id(int n)
{
  matrix<T> res(n,n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      res[i][j]= (i==j) ? 1 : 0;
  return res;
}

template<class T>
matrix<T> diag(int n, T first, ...)
{
  matrix<T> res=id<T>(n);
  T* ptr=&first;

  for(int i=0; i<n; i++)
  {
    res[i][i]=*(ptr);
    ptr++;
  }
  return res;
}

matrix<double> HilbertMatrix(int n);
matrix<double> RandomMatrix(int m, int n);

matrix<double> RandomSimilarMatrix(matrix<double> B);
column<double> RandomColumn(int n);
matrix<double> Sample1(int n);
matrix<double> Sample5(int n);
matrix<double> Sample8(int n);
matrix<double> zero(int m, int n);
column<double> zero(int n);

template<class T>
column<T> e(int n, int i)
{
  column<T> res(n);
  for(int k=0; k<n; k++) res[k]=1?i==k:0;
  return res;
}

// Returns matrices Linv, U that satisfy LU=A
template<class T>
void Gauss(matrix<T>& A, matrix<T>& Linv, matrix<T>& U)
{
  if(A.n!=A.m) throw string("non-square matrix");
  int n=A.n;
  matrix<T> P=A; // U
  matrix<T> C=id<T>(n); // Linv
  for(int j=0; j<n-1; j++)
  {
    int maxr=j;
    for(int r=j+1; r<n; r++)
    {
      maxr=(abs(P[maxr][j])<abs(P[r][j]))?r:maxr;
    }
    if(maxr!=j)
    {
      P.swap_row(j,maxr);
      C.swap_row(j,maxr);
    }
    for(int k=j+1; k<n; k++)
    {
      T t=-P[k][j]/P[j][j];
      P.plus_row(j, t, k);
      C.plus_row(j, t, k);
    }
  }
  Linv=C;
  U=P;
}

template<class T>
T sqr(T a)
{
  return a*a;
}

template<class T>
/**********************************************************\
  Solution x of the system Rx=b
  with R - right triangular matrix without zeros on diagonal
\**********************************************************/
void RTriangularSystem(matrix<T> R, column<T> b, column<T>& res_x)
{
  int N=R.n;
  column<T> x(N);
  for(int i=N-1; i>=0; i--)
  {
    T s=0;
    for(int j=i+1; j<N; j++) s+=R[i][j]*x[j];
    x[i]=(b[i]-s)/R[i][i];
  }
  res_x=x;
}

template<class T>
/**********************************************************\
  Solution x of the system Lx=b
  with L - left triangular matrix without zeros on diagonal
\**********************************************************/
void LTriangularSystem(matrix<T> L, column<T> b, column<T>& res_x)
{
  int N=L.n;
  column<T> x(N);
  for(int i=0; i<N; i++)
  {
    T s=0;
    for(int j=0; j<i; j++) s+=L[i][j]*x[j];
    x[i]=(b[i]-s)/L[i][i];
  }
  res_x=x;
}

// Returns matrices Qt,R that satisfy A=QR
template<class T>
void Householder(matrix<T>& A, matrix<T>& rQt, matrix<T>& rR)
/*
Householder's method
Qt*A=R
returns Qt:det Qt=1 and R
*/
{
  if(A.n!=A.m) throw string("non-square matrix");
  int n=A.n;
  int count=0;
  matrix<T> Qt=id<T>(n);
  matrix<T> R=A;
  column<T> u1(n), u2(n), u(n);
  row<T> r(n), p(n);
  for(int j=0; j<n-1; j++)
  {
    T normx=0;
    for(int i=j; i<n; i++) normx+=sqr(R[i][j]);
    normx=sqrt(normx);
    for(int i=j+1; i<n; i++) u1[i]=R[i][j];
    for(int i=j+1; i<n; i++) u2[i]=R[i][j];
    u1[j]=R[j][j]+normx; u2[j]=R[j][j]-normx;
    T normu1=0, normu2=0;
    for(int i=j; i<n; i++)
    {
      normu1+=sqr(u1[i]);
      normu2+=sqr(u2[i]);
    }
    for(int i=j; i<n; i++)
      if(normu1>normu2) u[i]=u1[i]; else u[i]=u2[i];
    T normu=sqrt(max(normu1,normu2));
    if(normu!=0)
    {
      for(int i=j; i<n; i++) u[i]/=normu;
      count++;
      // H*R
      for(int s=j; s<n; s++)
      {
        r[s]=0;
        for(int i=j; i<n; i++) r[s]+=u[i]*R[i][s];
      }
      for(int s=j; s<n; s++)
        for(int i=j; i<n; i++) R[i][s]-=2.0*u[i]*r[s];
      // H*Qt
      for(int s=0; s<n; s++)
      {
        p[s]=0;
        for(int i=j; i<n; i++) p[s]+=u[i]*Qt[i][s];
      }
      for(int s=0; s<n; s++)
        for(int i=j; i<n; i++) Qt[i][s]-=2.0*u[i]*p[s];
    }
  }
  /*if((count&1)==1)
  {
    for(int j=0; j<n; j++)
    {
      R[0][j]=-R[0][j];
      Qt[0][j]=-Qt[0][j];
    }
  }*/
  rR=R; rQt=Qt;
}

template<class T>
void Gradient_Descent(matrix<T>& A, column<T>& x, column<T>& b)
{
  column<T> r=A*x-b;
  double tau=(tp(r)*r)/(tp(A*r)*r);
  x=x-tau*r;
}

template<class T>
/***********************************
Input:
   A - square matrix with dominant eigenvalue
Output:
   res_lambda - dominant eigenvalue
   res_x - dominant eigenvector
***********************************/
void Power_Iteration(matrix<T>& A, double& res_lambda, column<T>& res_x, double eps)
{
  const int ITER_LIMIT=1000;
  column<T> x=RandomColumn(A.n);
  double lambda=0;
  int k=0;
  while (abs((tp(A*x)*x)/(tp(x)*x)-lambda)>=eps)
  {
    column<T> Ax=A*x;
    lambda=(tp(Ax)*x)/(tp(x)*x);
    x=Ax/Ax.norm();
    k++;
  }
  res_lambda=lambda;
  res_x=x;
  cout<<"Power Iteration. Number of iterations: "<<k<<endl;
}

template<class T>
/***********************************
Input:
   A - square matrix with dominant eigenvalue
Output:
   res_lambda - eigenvalue which is strictly nearest to zero
   res_x - dominant eigenvector
***********************************/
void Inverse_Power_Iteration(matrix<T>& A, double& res_lambda, column<T>& res_x, double eps)
{
  const int ITER_LIMIT=10000000;
  matrix<T> Qt,R;
  Householder(A, Qt, R);
  column<T> x=RandomColumn(A.n);
  x=x/x.norm();
  double lambda=0;
  int k=0;
  while (k<=ITER_LIMIT)
  {
    column<T> Aix;
    RTriangularSystem<double>(R, Qt*x, Aix);
    double next_lambda=(tp(x)*Aix)/(tp(Aix)*Aix);
    if(abs(lambda-next_lambda)<eps && k>0) break;
    x=Aix/Aix.norm();
    lambda=next_lambda;
    k++;
  }
  res_lambda=lambda;
  res_x=x;
  cout<<"Inverse Power Iteration. Number of iterations: "<<k<<endl;
}

template<class T>
/***********************************
Input:
   A - square matrix
Output:
   res_v - vector of eigenvalues
***********************************/
void QR_algorithm(matrix<T>& A, column<T>& res_v, matrix<T>& res_B, double eps)
{
  const int ITER_LIMIT=100000;
  matrix<T> A_cur, A_next;
  int N=A.n;
  column<T> v_cur(N), v_next(N);
  A_cur=A;  v_cur=A_cur.get_diag();
  int k=0;
  while (k<=ITER_LIMIT)
  {
    matrix<T> Qt,R;
    Householder(A_cur, Qt, R);
    A_next=R*tp(Qt); v_next=A_next.get_diag();
    if((v_next-v_cur).Cheb()<=eps) break;
    A_cur=A_next;  v_cur=A_cur.get_diag();
    k++;
  }
  res_B=A_cur;
  res_v=column<T>(N);
  for(int i=0; i<N; i++) res_v[i]=A_cur[i][i];
  cout<<"QR algorithm. Number of iterations: "<<k<<endl;
}

double energetic_norm(const matrix<double>& A, const column<double>& x);

column<double> getGridRegular(double a, double b, int N);



/********* STATIC VECTORS AND MATRICES *********
************************************************/

template<class T, int N>
class scolumn
{
public:
  T D[N];
  scolumn() { ;}
  ~scolumn() { ;}
  scolumn(const scolumn<T,N>& source)
  {
    for(int i=0; i<N; i++) D[i]=source[i];
  }
  T& operator[](int i)
  {
    return D[i];
  }
  const T& operator[](int i) const
  {
    return D[i];
  }
  void zero()
  {
    for(int i=0; i<N; i++)
    {
      D[i]=0;
    }
  }
  scolumn<T,N>& operator=(const scolumn<T,N>& S)
  {
    for(int i=0; i<N; i++)
    {
      D[i]=S[i];
    }
    return *this;
  }
  bool operator==(const scolumn<T,N>& S)
  {
    bool res=true;
    int j=0;
    while(res && j<N)
    {
      res=res && (D[j]==S[j]);
      j++;
    }
    return res;
  }
  bool operator!=(const scolumn<T,N>& S)
  {
    return !operator==(S);
  }
  scolumn<T,N> operator+(const scolumn<T,N>& S) const
  {
    scolumn<T,N> res;
    for(int i=0; i<N; i++)
    {
      res[i]=D[i]+S[i];
    }
    return res;
  }
  scolumn<T,N> operator-(const scolumn<T,N>& S) const
  {
    scolumn<T,N> res;
    for(int i=0; i<N; i++)
    {
      res[i]=D[i]-S[i];
    }
    return res;
  }
  scolumn<T,N> operator*(double a) const
  {
    scolumn<T,N> res;
    for(int i=0; i<N; i++)
    {
      res[i]=D[i]*a;
    }
    return res;
  }
  scolumn<T,N> operator/(double a) const
  {
    scolumn<T,N> res;
    for(int i=0; i<N; i++)
    {
      res[i]=D[i]/a;
    }
    return res;
  }
  T norm()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=D[i]*D[i];
    return sqrt(res);
  }
  T sqrnorm()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=D[i]*D[i];
    return res;
  }
  T norm1()
  {
    T res=0;
    for(int i=0; i<N; i++) res+=abs(D[i]);
    return res;
  }
  T Cheb()
  {
    T res=0;
	  for(int i=0; i<N; i++) res=max(res,abs(D[i]));
    return res;
  }
};

template<class T, int n>
scolumn<T,n> gets(const column<T>& S)
{
  scolumn<T,n> res;
  for(int i=0; i<n; i++)
  {
    res[i]=S[i];
  }
  return res;
}

template<class T, int m, int n>
class smatrix
{
  T M[m][n]; // col<row<T>>
public:
  smatrix()  {  ;}
  smatrix(const smatrix<T,m,n>& S)
  {
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      M[i][j]=S(i,j);
    }
  }
  ~smatrix()  { ;}
  T& operator()(int i, int j)
  {
    return M[i][j];
  }
  const T& operator()(int i, int j) const
  {
    return M[i][j];
  }
  void zero()
  {
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      M[i][j]=0;
    }
  }
  void id()
  {
    zero();
    for(int i=0; i<n; i++)
    {
      M[i][i]=1;
    }
  }
  smatrix<T,m,n>& operator=(const smatrix<T,m,n>& S)
  {
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      M[i][j]=S(i,j);
    }
    return *this;
  }
  smatrix<T,m,n> operator+(const smatrix<T,m,n>& S) const
  {
    smatrix<T,m,n> res;
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=M[i][j]+S(i,j);
    }
    return res;
  }
  smatrix<T,m,n> operator-(const smatrix<T,m,n>& S) const
  {
    smatrix<T,m,n> res;
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=M[i][j]-S(i,j);
    }
    return res;
  }
  smatrix<T,m,n> operator*(double a) const
  {
    smatrix<T,m,n> res;
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=M[i][j]*a;
    }
    return res;
  }
  smatrix<T,m,n> operator/(double a) const
  {
    smatrix<T,m,n> res;
    for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=M[i][j]/a;
    }
    return res;
  }
  scolumn<T,m> operator*(const scolumn<T,n>& x) const
  {
    scolumn<T,m> y;
    for(int i=0; i<m; i++)
    {
      y[i]=0;
      for(int j=0; j<m; j++)
      {
        y[i]+=M[i][j]*x[j];
      }
    }
    return y;
  }
  void swap_row(int i1, int i2)
  {
    for(int j=0; j<n; j++)
    {
      swap((*this)(i1,j), (*this)(i2,j));
    }
  }
  void plus_row(int i1, T ci1, int i2)
  {
    for(int j=0; j<n; j++)
    {
      (*this)(i2,j)+=(*this)(i1,j)*ci1;
    }
  }
};

template<class T, int m, int n>
smatrix<T,m,n> gets(const matrix<T>& S)
{
  smatrix<T,m,n> res;
  for(int i=0; i<m; i++)
  for(int j=0; j<n; j++)
  {
    res(i,j)=S(i,j);
  }
  return res;
}

template<class T, int m, int n>
scolumn<T,m> operator*(const smatrix<T,m,n>& A, const scolumn<T,n>& x)
{
  scolumn<T,n> res;
  for(int i=0; i<m; i++)
  {
    res[i]=0;
    for(int j=0; j<n; j++)
    {
      res[i]+=A(i,j)*x[j];
    }
  }
  return res;
}

// Warning! no size check
template<class T, int m, int n>
smatrix<T,m,n> operator*(const smatrix<T,m,n>& A, const smatrix<T,m,n>& B)
{
  smatrix<T,m,n> res;
  for(int i=0; i<m; i++)
  for(int j=0; j<n; j++)
  {
    res(i,j)=0;
    for(int k=0; k<n; k++)
    {
      res(i,j)+=A(i,k)*B(k,j);
    }
  }
  return res;
}

template<class T, int N>
scolumn<T,N> operator*(double a, const scolumn<T,N>& x)
{
  scolumn<T,N> res;
  for(int i=0; i<N; i++)
  {
    res[i]=a*x[i];
  }
  return res;
}

template<class T, int m, int n>
smatrix<T,m,n> operator*(double a, const smatrix<T,m,n>& A)
{
  smatrix<T,m,n> res;
  for(int i=0; i<m; i++)
  for(int j=0; j<n; j++)
  {
    res(i,j)=a*A(i,j);
  }
  return res;
}

// Warning! no size check
// Returns matrices Linv, U that satisfy LU=A
template<class T, int n>
void Gauss(const smatrix<T,n,n>& A, smatrix<T,n,n>& Linv, smatrix<T,n,n>& U)
{
  smatrix<T,n,n> P=A; // U
  smatrix<T,n,n> C; C.id(); // Linv
  for(int j=0; j<n-1; j++)
  {
    int maxr=j;
    for(int r=j+1; r<n; r++)
    {
      maxr=(abs(P(maxr,j))<abs(P(r,j)))?r:maxr;
    }
    if(maxr!=j)
    {
      P.swap_row(j,maxr);
      C.swap_row(j,maxr);
    }
    for(int k=j+1; k<n; k++)
    {
      T t=-P(k,j)/P(j,j);
      P.plus_row(j, t, k);
      C.plus_row(j, t, k);
    }
  }
  Linv=C;
  U=P;
}

template<class T, int n>
/**********************************************************\
  Solution x of the system Rx=b
  with R - right triangular matrix without zeros on diagonal
\**********************************************************/
void RTriangularSystem(smatrix<T,n,n> R, scolumn<T,n> b, scolumn<T,n>& res_x)
{
  scolumn<T,n> x;
  for(int i=n-1; i>=0; i--)
  {
    T s=0;
    for(int j=i+1; j<n; j++) s=s+R(i,j)*x[j];
    x[i]=(b[i]-s)/R(i,i);
  }
  res_x=x;
}

template<class T, int n>
smatrix<T,n,n> Inv(const smatrix<T,n,n>& A)
{
  smatrix<T,n,n> res,Linv,U; // LinvA=U
  Gauss(A,Linv,U);
  for(int k=0; k<n; k++)
  {
    scolumn<T,n> e;
    e.zero(); e[k]=1;
    scolumn<T,n> x;
    RTriangularSystem(U,Linv*e,x);
    for(int i=0; i<n; i++) res(i,k)=x[i];
  }
  return res;
}

template<class T, int n>
void SolveSLE(const smatrix<T,n,n>& A, const scolumn<T,n>& b, scolumn<T,n>& res_x)
{
  smatrix<T,n,n> Linv,R; // LinvA=U
  Gauss(A,Linv,R);
  RTriangularSystem(R,Linv*b,res_x);
}

template<class T, int n>
scolumn<T,n> SolveSLE(const smatrix<T,n,n>& A, const scolumn<T,n>& b)
{
  scolumn<T,n> res_x;
  smatrix<T,n,n> Linv,R; // LinvA=U
  Gauss(A,Linv,R);
  RTriangularSystem(R,Linv*b,res_x);
  return res_x;
}


template<class T, int N>
scolumn<T,N> zero()
{
  scolumn<T,N> res;
  for(int i=0; i<N; i++)
  {
    res[i]=0;
  }
  return res;
}

template<class T, int m, int n>
smatrix<T,m,n> zero()
{
  smatrix<T,m,n> res;
  for(int i=0; i<m; i++)
  for(int j=0; j<n; j++)
  {
    res(i,j)=0;
  }
  return res;
}


template<class T, int m, int n>
ostream& operator<<(ostream& out, smatrix<T,m,n> M)
{
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
    {
      out<<M(i,j)<<" ";
    }
    out<<endl;
  }
  return out;
}


template<class T, int n>
ostream& operator<<(ostream& out, scolumn<T,n> C)
{
  out<<"(";
  for(int i=0; i<n; i++)
  {
      out<<(i==0?"":",")<<C[i];
  }
  out<<")";
  return out;
}
template<int m, int n>
smatrix<double,m,n> RandomMatrix()
{
  //srand((unsigned int)time(NULL));
  smatrix<double,m,n> res;
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=(rand()%2000)/1000.0-1.0;
    }
  return res;
}

template<class T, int m, int n>
smatrix<T,m,n> tp(const smatrix<T,m,n>& S)
{
  smatrix<T,m,n> res;
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=S(j,i);
    }
  return res;
}

template<class T, int m, int n>
smatrix<T,m,n> trunc(const smatrix<T,m,n>& S)
{
  smatrix<T,m,n> res;
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
    {
      res(i,j)=int(S(i,j));
    }
  return res;
}


template<class T>
std::vector<T> get_vector(const column<T>& A)
{
  std::vector<T> res;
  for(int i=0; i<A.N; i++)
  {
    res.push_back(A[i]);
  }
  return res;
}

template<class T, int N>
bool check_validity(const scolumn<T,N>& V)
{
  bool res=true;
  for(int i=0; i<N; i++)
  {
    res = res && check_validity(V[i]);
  }
  return res;
}
