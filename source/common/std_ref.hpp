#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cassert>

using std::string;
using std::istream;
using std::ifstream;
using std::ostream;
using std::ofstream;
using std::ostringstream;
using std::istringstream;
using std::max;
using std::cout;
using std::cin;
using std::min;
using std::swap;
using std::endl;
using std::stringstream;
using std::vector;
using std::abs;

template<class T>
string convert_to_string(const T& U)
{
  stringstream buf;
  buf.str("");
  buf<<U;
  string res;
  buf>>res;
  return res;
}
