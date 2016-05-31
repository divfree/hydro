#pragma once

#include "../common/std_ref.hpp"
#include "../common/data_structures.hpp"
#include <boost/thread/thread.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include "experiment.hpp"
#include "pvar.hpp"

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;

#define CMap_number 4
extern string Map_name[Map_number];

void sleep(int seconds);

class TConsole
{
public:
  typedef void (TConsole::*TCmd)(string);
  TConsole();
  ~TConsole();
  BinSearchSet<string, TCmd> Commands;
  BinSearchSet<string, string> Commands_desc;
  void init_commands(TConsole& console);
  void error(string msg);
  void write_msg(string msg);
  void msg_eol(string msg);
  string prompt();
  string cur_exp_name;
  // COMMANDS
  void cmd_exit(string arg);
  void cmd_echo(string arg);
  void cmd_add_experiment(string arg);
  void cmd_current_experiment(string arg);
  void cmd_del_experiment(string arg);
  void cmd_init(string arg);
  void cmd_start(string arg);
  void cmd_start_in_main_thread(string arg);
  void cmd_stop(string arg);
  void cmd_list_exp(string arg);
  void cmd_set(string arg);
  void cmd_list_parameters(string arg);
  void cmd_del(string arg);
  void cmd_test(string arg);
  void cmd_test2(string arg);
  void cmd_help(string arg);
  void cmd_run(string arg);
  void cmd_sleep(string arg);
  void cmd_usleep(string arg);
  void cmd_term(string arg);
  void cmd_nop(string arg);
  void cmd_cd(string arg);
  void cmd_dir(string arg);
  void cmd_save(string arg);
  void cmd_load(string arg);
  //void cmd_more_time(string arg);
  void cmd_draw_mu_mix(string arg);
  void cmd_hello(string arg);
  void cmd_mkdir(string arg);
  void cmd_wait_for_experiments_completion(string arg);
  void cmd_save_grid_double(string arg);
  void cmd_load_grid_double(string arg);
  void cmd_copy_U_comp_to_grid_double(string arg);
  void cmd_status(string arg);
  void cmd_pvar_new(string arg);
  void cmd_pvar_start(string arg);
  void cmd_pvar_set(string arg);
  void cmd_pvar_table(string arg);
  void cmd_value(string arg);
  void cmd_open_log(string arg);
  void cmd_write_parabolic_profile_plt(string arg);
  void cmd_write_pois_cyl(string arg);

  ofstream cflog;

  // EXP CONTROL FUNCTIONS
  void exp_start(TExperiment*);
  void exp_create_thread(TExperiment*);
  void exp_del(TExperiment*);
  void exp_term(TExperiment*);
  void exp_stop(TExperiment*);

  TExperiment* check_cur_exp();
  void execute();
  bool flag_exit;
  bool finished;
  BinSearchSet2<string, bool> CP_bool;
  BinSearchSet2<string, double> CP_double;
  BinSearchSet2<string, int> CP_int;
  BinSearchSet2<string, string> CP_string;
  TParameter<string>* CMap[CMap_number];
  BinSearchSet<string, TExperiment*> Experiments;
  void open_file(ofstream& fout, string filename);
  void open_file(ifstream& fin, string filename);
  void check_directory(string dirname);
  void check_file_exists(string filename);
  void check_file_good(ofstream& file, string filename);
  void check_file_good(ifstream& file, string filename);
  int running_exp_count();
  int undone_exp_count();
  int available_threads_number();
  void vary_grid_check(ofstream& vg_out);
  void vary_grid_step_check(ofstream& vg_out);
  void confirmation_yes_no(string msg);

  void scheduler_init();
  void scheduler_thread();
  boost::thread* scheduler_thread_ptr;
  void scheduler_term();
  bool scheduler_terminate;
  bool scheduler_terminate_done;
  int threads_count;
  int pending_count;

  PVar* pvar;

  string get_parameter_value(string cmdline, TExperiment* ex=0);
  string pname_parser(string cmdline);
  bool valid_char_for_parameter_name(char c);
};

void extract_cmd_name_and_arg(string str, string& name, string& arg);
