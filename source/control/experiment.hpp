#pragma once

#include "../common/std_ref.hpp"
#include "parameters.hpp"
#include "../common/linear_algebra.hpp"
#include "../common/data_structures.hpp"
#include "../common/common.hpp"
#include <thread>
#include <chrono>
#include "../control/metrics.hpp"
#include <memory>
#include <map>

#define Map_number 5
extern string Map_name[Map_number];

class TConsole;
class TModule;
class ModuleFactory;

namespace hydro2D_uniform
{
  class hydro;
}

namespace hydro3D_structured
{
  class hydro;
}

enum class ES { blank, progress_init, done_init, pending, thread_created, progress_step, done_step, progress_write_results, done_write_results, thread_destroyed, finished, destroyed, error, terminated };

class TExperiment
{
public:
  TExperiment(TConsole* _console);
  ~TExperiment();
  TConsole* console;
  ES status;
  void status_change(ES new_status);
  string status_name(ES st);
  string status_name();
  string name;
  string output_dir;
  bool st_thread, st_init, st_term, st_pending, st_error, st_finished;
  bool sig_term;

  // P Parameters
  BinSearchSet2<string, int> P_int;
  BinSearchSet2<string, double> P_double;
  BinSearchSet2<string, string> P_string;
  BinSearchSet2<string, bool> P_bool;
  BinSearchSet2<string, column<double>> P_vect;
  TParameter<string>* Map[Map_number];
  bool flag(string flag_name);
  bool uflag(string flag_name); // flag() without existence check, direct access by name

  void thread();
  std::shared_ptr<std::thread> thread_ptr;

  void open_file(ofstream& fout, string filename);
  void open_file(ifstream& fin, string filename);
  ofstream flog;
  void logger(string msg);

  std::shared_ptr<TModule> module;

  MultiTimer<std::string> timer_;

  // CONTROL FUNCTIONS
  void init();
  void term();

private:
  static std::map<std::string, std::shared_ptr<ModuleFactory>> module_map_;
  enum class MODULE {hydro2D_uniform, hydro3D_structured, proppant2D};
  void set_name();
  // Subsystems init/term functions
  void init_log();

public:
  static void RegisterModule(std::string name,
                             std::shared_ptr<ModuleFactory> factory) {
    module_map_[name] = factory;
  }
};

class TExperiment_ref
{
  public:
  TExperiment* ex;
  TExperiment_ref(TExperiment* _ex) : ex(_ex),
  P_bool(ex->P_bool), P_int(ex->P_int), P_double(ex->P_double), P_string(ex->P_string), P_vect(ex->P_vect), output_dir(ex->output_dir), flog(ex->flog) {;}
  BinSearchSet2<string, bool>& P_bool;
  BinSearchSet2<string, int>& P_int;
  BinSearchSet2<string, double>& P_double;
  BinSearchSet2<string, string>& P_string;
  BinSearchSet2<string, column<double>>& P_vect;
  inline bool flag(string flag_name) { return ex->flag(flag_name); }
  inline void status_change(ES new_status) { ex->status_change(new_status); }
  string& output_dir;
  ofstream& flog;
};

class TModule : virtual public TExperiment_ref
{
  protected:
  void increase_time();
  double get_increased_time(double t, double dt);

  public:
  TModule(TExperiment* _ex);
  virtual ~TModule() {;}
  void thread();
  virtual void write_results(bool force=false)=0;
  void update_dt();
  double dt;
  double tn;
  double tnp;
  virtual void save_state(std::ofstream&) { cout<<"Warning: save_state() function is not defined for this module"<<endl;}
  virtual void load_state(std::ifstream&) { cout<<"Warning: load_state() function is not defined for this module"<<endl;}
  void save_state(string filename);
  void load_state(string filename);

  private:
  virtual void step()=0;
  void cycle();
  void init_step_time();
  void update_step_time(double seconds);
  void write_step_header(double t);
  void write_step_header();
  void write_step_footer();
};

class ModuleFactory {
 public:
  virtual std::shared_ptr<TModule> Create(TExperiment* ex) const = 0;
  virtual ~ModuleFactory() {}
};

template <class Module>
class ModuleFactoryTemplate : public ModuleFactory {
 public:
  std::shared_ptr<TModule> Create(TExperiment* ex) const override {
    return std::make_shared<Module>(ex);
  }
};

template <class Module>
class ModuleRegistrator {
 public:
  ModuleRegistrator(std::vector<std::string> name_aliases) {
    std::shared_ptr<ModuleFactory> factory(
        new ModuleFactoryTemplate<Module>());
    for (auto name : name_aliases) {
      TExperiment::RegisterModule(name, factory);
    }
  }
};

class Iterations : virtual private TExperiment_ref
{
  protected:
  Iterations(TExperiment* _ex) : TExperiment_ref(_ex) {;}
  void iter_history_open();
  void iter_history_write(int s);
  bool iter_history_condition();
  void write_stat_s(int s, double R);
  void iter_history_close();
  bool while_condition(int s, double R);
  virtual void iter_history_write_scalar(ofstream& fout)=0;
  virtual void iter_history_variables(ofstream& fout)=0;
  virtual void iter_history_mesh_write(ofstream&) {;}
  virtual void iter_history_mesh_header(ofstream&) {;}
  private:
  ofstream fihxy, fihs;
};

// USAGE: corresponding handlers will be called on certain events
// (hard implementation, i.e. handlers are embedded in the class)
class signal_analyzer : private TExperiment_ref
{
public:
  int count;
  signal_analyzer(TExperiment* _ex, double _relax=1.0, double _min_diff_threshold=1E-3);
  double t, t_prev;
  double uc, uc_prev; // u corrected
  double relax;  // 0..1
  void step(double t, double u);
  double min_value, min_interval, min_t, min_value_diff_rel, min_diff_threshold;
  bool descending;
  void min_handler();
  void min_diff_threshold_handler();
  bool threshold_reached;
};
