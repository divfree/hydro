#include "experiment.hpp"

#include "../control/console.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "../revision.hpp"


std::map<std::string, std::shared_ptr<ModuleFactory>> TExperiment::module_map_;

string Map_name[Map_number]={"int", "double", "string", "bool", "vect"};

void TExperiment::status_change(ES new_status)
{
  status=new_status;
  if(status==ES::error)
  {
    cout<<"Status of '"<<name<<"' ("<<ecast(P_string("name"))<<") changed to '"<<status_name(status)<<"'"<<endl;
  }
  #ifdef _DEBUG2
  cout<<"STATUS: "<<status_name()<<endl;
  #endif
}

string TExperiment::status_name(ES status) {
  switch(status) {
    case ES::blank:
      return "blank";
    case ES::progress_init:
      return "progress_init";
    case ES::done_init:
      return "done_init";
    case ES::pending:
      return "pending";
    case ES::thread_created:
      return "thread_created";
    case ES::progress_step:
      return "progress_step";
    case ES::done_step:
      return "done_step";
    case ES::progress_write_results:
      return "progress_write_results";
    case ES::done_write_results:
      return "done_write_results";
    case ES::finished:
      return "finished";
    case ES::thread_destroyed:
      return "thread_destroyed";
    case ES::error:
      return "error";
    default:
      return "Unknown_status";
  }
}

string TExperiment::status_name() {
  return status_name(status);
}

TExperiment::TExperiment(TConsole* _console) : console(_console) {
  thread_ptr=0;

  TParameter<string>* Map1[Map_number]={&P_int, &P_double, &P_string, &P_bool, &P_vect};
  for(int k=0; k<Map_number; k++)
  {
    Map[k]=Map1[k];
  }

  st_init=false;
  st_term=false;
  st_thread=false;
  st_pending=false;
  st_error=false;
  st_finished=false;
  sig_term=false;

  status_change(ES::blank);
  P_string.set("status", status_name());

#ifdef _OPENMP
  P_int.set("omp_num_threads", omp_get_max_threads());
#else
  P_int.set("omp_num_threads", 1);
#endif

  set_name();
}

void TExperiment::init() {
  output_dir=ecast(P_string(_output_dir));
  if(output_dir!="" && !directory_exists(output_dir))
  {
    create_directory(output_dir);
  }

  init_log();

  // MODULE selection
  string module_name=P_string.exist("MODULE") ? P_string["MODULE"]
                                              : "hydro2D_uniform_MPI";

  if (module_map_.count(module_name)) {
    module = module_map_[module_name]->Create(this);
  } else {
    throw string("Unknown module '"+module_name+"'");
  }

  st_init=true;
}

void TExperiment::set_name()
{
  P_string.set(_exp_name, ""+ecast(P_string(_exp_name)));
}

bool TExperiment::flag(string flag_name)
{
  bool* ptr=P_bool(flag_name);
  return (ptr && *ptr) ;
}

bool TExperiment::uflag(string flag_name)
{
  return P_bool[flag_name];
}

void TExperiment::open_file(ofstream& fout, string filename)
{
  console->open_file(fout, output_dir+filename);
}

void TExperiment::open_file(ifstream& fin, string filename)
{
  console->open_file(fin, output_dir+filename);
}

void TExperiment::init_log()
{
  if(!ecast(P_bool("no_output")))
  {
    string filename_log=*P_string.set("filename_log", P_string[_exp_name]+".log");
    open_file(flog, filename_log);
    flog << "# hydro3 experiment\n";
    flog << "# revision " << kGitRevision << "\n";
    flog << "# EXPERIMENT PARAMETERS" << std::endl;

    for(int k=0; k<Map_number; k++)
    if(Map[k]->get_length()>0)
    {
      flog<<"# "<<Map_name[k]<<" parameters"<<endl;
      for(int i=0; i<Map[k]->get_length(); i++)
      {
        flog<<"set "<<Map_name[k]<<" "<<Map[k]->get_key(i)<<" ";
        Map[k]->write_data_by_index(flog, i);
        flog<<endl;
      }
      flog<<endl;
    }
  }
}

void TExperiment::term()
{
  string* after_finish_ptr=P_string("after_finish_script");
  if(after_finish_ptr)
  {
    string script=*after_finish_ptr;
    console->cmd_run(script);
  }

  if(!ecast(P_bool("no_output")))
  {
    flog<<endl<<"Experiment terminated."<<endl;
    double /*t_FVM=P_double["total_time_FVM"], t_part=P_double["total_time_particles"], */t_all=P_double["total_step_time"];
    flog<<"Execution time:"<<endl;
    /*flog<<"FVM: "<<seconds_to_hms(t_FVM)<<" ("<<round(t_FVM/t_all*10000)/100.0<<"%)"<<endl;
    flog<<"particles: "<<seconds_to_hms(t_part)<<" ("<<round(t_part/t_all*10000)/100.0<<"%)"<<endl;*/
    flog<<"all: "<<seconds_to_hms(t_all)<<endl;

    if(P_int.exist("cells_number") && P_int.exist("s_sum"))
    {
      P_double.set("t_cell", t_all/(P_int["cells_number"]*P_int["s_sum"]));
      flog<<"one cell, one iter: "<<P_double["t_cell"]<<endl;
    }

    flog << std::endl;
    for (auto entry : timer_.GetTotalTime()) {
      flog << entry.first << " : " << entry.second << std::endl;
    }

    flog.close();
  }

  console->pvar->table_record(this);

  st_term=true;
}

void TExperiment::thread()
{
  module->thread();
}

void TExperiment::logger(string msg)
{
  flog<<msg<<endl;
}

TExperiment::~TExperiment()
{

}
