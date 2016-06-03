/*******************************************************************
*****************    Parameters Variator (PVar)     ****************
********************************************************************/

#pragma once

#include "../common/common.hpp"
#include "../common/linear_algebra.hpp"
#include <thread>
#include <memory>
#include <atomic>

class TConsole;
class TExperiment;

class PVar_single
{
public:
  PVar_single();
  PVar_single(string arg);
  string pname;
  string ptype; // parameter type
  string vtype; // variator type
  column<string> values;
  string current;
  void nv_step_add(stringstream& buf);
  void nv_step_mul(stringstream& buf);
  void nv_range_add(stringstream& buf);
  void nv_range_mul(stringstream& buf);
  void nv_enum(stringstream& buf);
};

class PVar
{
public:
  PVar(TConsole* parent);
  TConsole* console;
  column<PVar_single> VL; // variators list
  void new_variator(string arg);
  void init();
  void start();
  void thread();
  void step(int depth);
  void term();
  std::shared_ptr<std::thread> thread_ptr;
  int variations_count;
  std::atomic<bool> running;
  ofstream ftable;
  column<string> table_columns;
  void table_open(string arg);
  void table_record(TExperiment* ex=0);
  void table_close();
};

#include "console.hpp"
#include "experiment.hpp"
