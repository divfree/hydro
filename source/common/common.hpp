#pragma once

#include <time.h>
#include <stdio.h>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/round.hpp>
#include "std_ref.hpp"

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;

string get_time(string fmt);
string get_string_with_ws(istream& in);
string get_current_directory();
void set_current_directory(string path);
string IntToStr(int a);
string DoubleToStr(double a, int precision=5);
double StrToDouble(string s);
bool directory_exists(string dirname);
void create_directory(string dirname);
void remove_file(string filename);
bool file_exists(string filename);
string point_to_string(int i, int j);

string ecast(string* ptr);
int ecast(int* ptr, int def=0);
bool ecast(bool* ptr);
double ecast(double* ptr, double def=0.0);
int round05(double a);
//int round(double a);
int myround(double a);

int simple_pow(int x, int y);
double sqr(double a);

string seconds_to_hms(double sec);

// return k : b=10^k-1, b>=a, k->min
int sup_power10_k(int a);
int sup_power10_k(double a);

double round_to_exponent(double a, int k);

double round_to_precision(double a, int precision);

double max(double, double, double);

struct TRect
{
  int ai, aj, bi, bj;
  TRect() { }
  TRect(int _ai, int _aj, int _bi, int _bj)
  {
    ai=_ai;
    aj=_aj;
    bi=_bi;
    bj=_bj;
  }
  TRect xm()
  {
    return TRect(ai,aj,bi+1,bj);
  }
  TRect ym()
  {
    return TRect(ai,aj,bi,bj+1);
  }
  TRect expand_by_1()
  {
    return TRect(ai-1,aj-1,bi+1,bj+1);
  }
};

bool in_rect(int i, int j, int i1, int j1, int i2, int j2);
bool in_rect(int i, int j, const TRect& rect);

namespace sysinfo
{
  int virtual_usage_kb();
  int physical_usage_kb();
  int threads();
}

double get_seconds_count(cpu_timer& timer);

double pow3(double a);
double pow4(double a);
double bilinear_interpolation(double x, double y, double x00, double y00, double x11, double y11, double a00, double a10, double a11, double a01);
double linear_interpolation(double x, double x0, double x1, double a0, double a1);
double distance2(double x1, double y1, double x2, double y2);

template <class TData>
void read_data(istream& in, TData& data)
{
  in>>data;
}

void read_data(istream& in, string& data);
