#pragma once

#include "linear_algebra.hpp"
#include "../common/common.hpp"

constexpr double pi=atan(1.0)*4;

namespace hydro2D_uniform
{

typedef column<double> vect;
typedef column<int> arrayint;
typedef matrix<double> matr;

typedef grid<double> grid_double;

#define U_comp_u 0
#define U_comp_v 1
#define U_comp_p 2
#define U_comp_beta 3

// static vector experiment
struct svectexp : public scolumn<double, 4>
{
  svectexp() : u(D[U_comp_u]), v(D[U_comp_v]), p(D[U_comp_p]), beta(D[U_comp_beta]) {}
  void operator=(const svectexp& V)
  {
    scolumn::operator=(V);
  }
  svectexp(const scolumn& V) : u(D[U_comp_u]), v(D[U_comp_v]), p(D[U_comp_p]), beta(D[U_comp_beta])
  {
    scolumn::operator=(V);
  }
  svectexp(const column<double>& V) : u(D[U_comp_u]), v(D[U_comp_v]), p(D[U_comp_p]), beta(D[U_comp_beta])
  {
    for(int i=0; i<4; i++)
    {
      D[i]=V[i];
    }
  }
  double& u, &v, &p, &beta;
};

// static vector spatial
struct svectsp : public scolumn<double, 2>
{
  svectsp() : x(D[0]), y(D[1]) {}
  void operator=(const svectsp& V)
  {
    scolumn::operator=(V);
  }
  svectsp(const scolumn& V) : x(D[0]), y(D[1])
  {
    scolumn::operator=(V);
  }
  svectsp(const column<double>& V) : x(D[0]), y(D[1])
  {
    for(int i=0; i<2; i++)
    {
      D[i]=V[i];
    }
  }
  double& x, &y;
};

typedef grid<svectexp> grid_svectexp;

matr get_matr(double a11, double a12, double a13, double a14, double a21, double a22, double a23, double a24,
              double a31, double a32, double a33, double a34, double a41, double a42, double a43, double a44);
matr get_matr(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);
vect get_vect(double a1, double a2, double a3, double a4);
vect get_vect(double a1, double a2, double a3);
vect get_vect(double a1, double a2);

double norm_L1_mean(const grid_double& U, double hx, double hy);
double distance2(double x1, double y1, double x2, double y2);
//void init_array(grid_svectexp& U, int Nx, int Ny, svectexp value);
//void init_array(grid_double& U, int Nx, int Ny, double value);

template<class T>
void init_array(grid<T>& U, int Nx, int Ny, T value)
{
  U.change_N(Nx,Ny);
  for(int i=0; i<Nx; i++)
  {
    for(int j=0; j<Ny; j++)
    {
      U(i,j)=value;
    }
  }
}

scolumn<double,3> get_svect(double a1, double a2, double a3);
smatrix<double,3,3> get_smatr(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);

scolumn<double,1> get_svect(double a1);
smatrix<double,1,1> get_smatr(double a11);

double get_double(smatrix<double,1,1> U);
double get_double(scolumn<double,1> U);

template <int sdim>
void write_comp(const grid<scolumn<double,sdim>>& U, int k)
{
  for(int j=U.Ny; j>=0; j--)
  {
    for(int i=0; i<U.Nx; i++)
    {
      cout<<U(i,j)[k]<<" ";
    }
    cout<<endl;
  }
}

template <int sdim>
double norm_max(const grid<scolumn<double,sdim>>& U, int k)
{
  double res=0;
  int Nx=U.Nx, Ny=U.Ny;
  for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++)
  {
    res=max(res, abs(U(i,j)[k]));
  }
  return res;
}

template <int sdim>
double norm_L1_mean(const grid<scolumn<double,sdim>>& U, int k, double hx, double hy)
{
  double res=0;
  for(int i=0; i<U.Nx; i++)
  for(int j=0; j<U.Ny; j++)
  {
    res+=abs(U(i,j)[k])*hx*hy;
  }
  return res/(U.Nx*hx*U.Ny*hy);
}

template <int sdim>
void get_vect_component(grid_double& u, grid<scolumn<double,sdim>>& Q, int k)
{
  int Nx=Q.Nx;
  int Ny=Q.Ny;
  init_array(u, Nx, Ny, 0.0);

  for(int i=0; i<Nx; i++)
  {
    for(int j=0; j<Ny; j++)
    {
      u(i,j)=Q(i,j)[k];
    }
  }
}

void grid_double_to_grid_svectexp(const grid_double& u, grid_svectexp& U, int component);
void grid_svectexp_to_grid_double(grid_double& u, const grid_svectexp& U, int component);


// grid storing format: Nx, Ny, EOF
// then Ny lines of DATA (j=Ny-1..0, i=0..Nx-1)
template<class T>
void save_grid(const grid<T>& u, ostream& out)
{
  out<<u.Nx<<" "<<u.Ny<<endl;
  for(int j=u.Ny-1; j>=0; j--)
  {
    for(int i=0; i<u.Nx; i++)
    {
      out<<u(i,j)<<" ";
    }
    out<<endl;
  }
}

template<class T>
void load_grid(grid<T>& u, istream& in)
{
  int iNx, iNy;
  in>>iNx>>iNy;
  u.change_N(iNx,iNy);
  for(int j=u.Ny-1; j>=0; j--)
  {
    for(int i=0; i<u.Nx; i++)
    {
      double a=get_nan();
      in>>a;
      if(!check_validity(a)) throw string("load_grid: invalid or missing number for "+point_to_string(i,j)+" node");
      u(i,j)=a;
    }
  }
}

template<class T>
void write(const grid<T>& u, string name)
{
  cout<<"GRID: "<<name<<endl<<u<<endl;
}

// namespace end
}
