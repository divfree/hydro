#include "linear_algebra_ext.hpp"

namespace hydro2D_uniform
{

matr get_matr(double a11, double a12, double a13, double a14,
               double a21, double a22, double a23, double a24,
               double a31, double a32, double a33, double a34,
               double a41, double a42, double a43, double a44)
{
  matr A(4,4);
  A[0][0]=a11; A[0][1]=a12; A[0][2]=a13; A[0][3]=a14;
  A[1][0]=a21; A[1][1]=a22; A[1][2]=a23; A[1][3]=a24;
  A[2][0]=a31; A[2][1]=a32; A[2][2]=a33; A[2][3]=a34;
  A[3][0]=a41; A[3][1]=a42; A[3][2]=a43; A[3][3]=a44;
  return A;
}

matr get_matr(double a11, double a12, double a13,
               double a21, double a22, double a23,
               double a31, double a32, double a33)
{
  matr A(3,3);
  A[0][0]=a11; A[0][1]=a12; A[0][2]=a13;
  A[1][0]=a21; A[1][1]=a22; A[1][2]=a23;
  A[2][0]=a31; A[2][1]=a32; A[2][2]=a33;
  return A;
}

vect get_vect(double a1, double a2, double a3, double a4)
{
  vect a(4);
  a[0]=a1; a[1]=a2; a[2]=a3; a[3]=a4;
  return a;
}

vect get_vect(double a1, double a2, double a3)
{
  vect a(3);
  a[0]=a1; a[1]=a2; a[2]=a3;
  return a;
}

vect get_vect(double a1, double a2)
{
  vect a(2);
  a[0]=a1; a[1]=a2;
  return a;
}

smatrix<double,3,3> get_smatr(double a11, double a12, double a13,
               double a21, double a22, double a23,
               double a31, double a32, double a33)
{
  smatrix<double,3,3> A;
  A(0,0)=a11; A(0,1)=a12; A(0,2)=a13;
  A(1,0)=a21; A(1,1)=a22; A(1,2)=a23;
  A(2,0)=a31; A(2,1)=a32; A(2,2)=a33;
  return A;
}

scolumn<double,3> get_svect(double a1, double a2, double a3)
{
  scolumn<double,3> a;
  a[0]=a1; a[1]=a2; a[2]=a3;
  return a;
}

smatrix<double,1,1> get_smatr(double a11)
{
  smatrix<double,1,1> A;
  A(0,0)=a11;
  return A;
}

scolumn<double,1> get_svect(double a1)
{
  scolumn<double,1> a;
  a[0]=a1;
  return a;
}

double get_double(smatrix<double,1,1> U)
{
  return U(0,0);
}

double get_double(scolumn<double,1> U)
{
  return U[0];
}

double norm_max(const grid_double& U)
{
  double res=0;
  int Nx=U.Nx, Ny=U.Ny;
  for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++)
  {
    res=max(res, abs(U(i,j)));
  }
  return res;
}

double norm_L1_mean(const grid_double& U, double hx, double hy)
{
  double res=0;
  for(int i=0; i<U.Nx; i++)
  for(int j=0; j<U.Ny; j++)
  {
    res+=abs(U(i,j))*hx*hy;
  }
  return res/(U.Nx*hx*U.Ny*hy);
}

void grid_double_to_grid_svectexp(const grid_double& u, grid_svectexp& U, int k)
{
  int Nx=u.Nx, Ny=u.Ny;
  U.change_N(Nx,Ny);
  for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++)
  {
    U(i,j)[k]=u(i,j);
  }
}

void grid_svectexp_to_grid_double(grid_double& u, const grid_svectexp& U, int k)
{
  int Nx=U.Nx, Ny=U.Ny;
  u.change_N(Nx,Ny);
  for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++)
  {
    u(i,j)=U(i,j)[k];
  }
}

// namespace end
}
