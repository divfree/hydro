#include "experiment.hpp"
#include "console.hpp"
#include "metrics.hpp"

void TModule::init_step_time()
{
  ex->P_double.set("total_step_time", 0);
  ex->P_double.set("last_step_time", 0);
  ex->P_double.set("memory", 0);
}

TModule::TModule(TExperiment* _ex) : TExperiment_ref(_ex)
{
  dt=0.0;
  tn=0.0;
  P_int.set("n", 0);
  P_double.set("t", tn);
  update_dt();
  tnp=get_increased_time(tn,dt);
  init_step_time();
}

void TModule::update_step_time(double seconds) {
  ex->P_double["last_step_time"] = seconds;
  ex->P_double["total_step_time"] += seconds;
  ex->P_double["memory"] = sysinfo::physical_usage_kb()/1000.0;
}


void TModule::update_dt()
{
  double dt_old=dt;
  if(P_double.exist("T0") && P_double.exist("dt0") && (tn<P_double["T0"])) {
    dt=P_double["dt0"];
  }
  else if(P_double.exist("T1") && P_double.exist("dt1") && (tn<P_double["T1"]))  {
    dt=P_double["dt1"];
  }
  else dt=P_double["dt"];

  if(dt!=dt_old) flog<<"Time step changed from "<<dt_old<<" to "<<dt<<endl;
}

void TModule::cycle()
{
  if(!ex->st_finished)
  {
    SingleTimer timer;

    write_step_header();

    step();
    increase_time();

    if(tn>=P_double["T"])
    {
      ex->st_finished=true;
    }

    write_results(ex->st_finished);

    update_step_time(timer.GetSeconds());

    write_step_footer();

    if(P_double.exist("eps_n") && P_double.exist("last_Rn") && P_double["last_Rn"]<P_double["eps_n"])
    {
      write_results(true);
      P_double["T"]=tn;
      ex->st_finished=true;
    }
  }
}

void TModule::increase_time()
{
  P_int["n"]++;
  tn=tnp;
  P_double["t"]=tn;
  update_dt();
  tnp=get_increased_time(tn,dt);
}

double TModule::get_increased_time(double t, double dt)
{
  //return round_to_exponent(t+dt, sup_power10_k(1.0/dt)+2);
  return t+dt;
}

void TModule::thread()
{
  try
  {
    while(!ex->st_finished && !ex->sig_term)
    {
      cycle();
    }
    ex->term();
    cout<<"("<<ex->name<<", "<<ex->P_string["name"]<<" ) done with status "<<ex->status_name()<<endl;
    ex->console->cflog<<"("<<ex->name<<", "<<ex->P_string["name"]<<" ) done with status "<<ex->status_name()<<", "<<ex->console->undone_exp_count()-1<<" left"<<endl;
  }
  catch(string msg)
  {
    ex->st_error=true;
    cout<<"ERROR: "<<msg<<endl;
    ex->flog<<"ERROR: "<<msg<<endl;
  }
  ex->st_thread=false;
  ex->console->threads_count--;
}


void TModule::write_step_header()
{
  ex->flog<<"BEGIN. t="<<tnp<<", n="<<P_int["n"]+1<<endl;
}

void TModule::write_step_footer()
{
  ex->flog<<"END. stat: s="<<P_int["last_s"]<<", Rs="<<P_double["last_R"];
  if(P_double.exist("last_Rn")) ex->flog<<", Rn="<<P_double["last_Rn"];
  ex->flog<<", t_all="<<P_double["last_step_time"]<<", mem="<<P_double["memory"]<<"MB"<<endl<<endl;
}

void TModule::save_state(string filename)
{
  ofstream fsave;
  ex->open_file(fsave, filename);
  fsave.precision(16);
  save_state(fsave);
  fsave.close();
}

void TModule::load_state(string filename)
{
  ifstream fload;
  ex->open_file(fload, filename);
  load_state(fload);
  fload.close();
}

void Iterations::write_stat_s(int s, double R)
{
  if(ex->flag(_stat_s_enable))
  {
    flog<<".....s="<<s<<", Rs="<<R<<endl;
  }
}

void Iterations::iter_history_open()
{
  if(iter_history_condition())
  {
    ex->open_file(fihs, P_string[_exp_name]+"_iter_history.scalar.plt");
    fihs<<"VARIABLES=";
    iter_history_variables(fihs);
    fihs<<endl;
    fihs<<"ZONE ZONETYPE=ORDERED DATAPACKING=POINT"<<endl;
    fihs<<"T=\""<<P_string[_plt_title]<<"\""<<endl;
    fihs<<"# History of iterations"<<endl;

    if(ex->flag(_iter_history_mesh))
    {
      ex->open_file(fihxy, P_string[_exp_name]+"_iter_history.grid.plt");
      iter_history_mesh_header(fihxy);
      fihxy<<endl;
    }
  }
}

void Iterations::iter_history_write(int s)
{
  if(iter_history_condition())
  {
    P_int["s"]=s;

    iter_history_write_scalar(fihs);
    fihs<<endl;

    if(ex->flag(_iter_history_mesh) && s%P_int["iter_history_mesh_skip"]==0)
    {
      iter_history_mesh_write(fihxy);
      fihxy<<endl;
    }
  }
}

bool Iterations::iter_history_condition()
{
  return (P_bool["iter_history_enable"] && P_int["iter_history_n"]==P_int["n"]+1);
}

void Iterations::iter_history_close()
{
  if(iter_history_condition())
  {
    fihs.close();

    if(ex->flag(_iter_history_mesh))
    {
      fihxy.close();
    }
  }
}

bool Iterations::while_condition(int s, double R)
{
  // iter_history enabled
  if(iter_history_condition() && P_int("iter_history_sfixed") && P_int["iter_history_sfixed"]>0)
  {
    if(s>=P_int["iter_history_sfixed"]) return false; else return true;
  }
  else
  // sfixed enabled
  if(P_int["sfixed"]>=0)
  {
    if(s>=P_int["sfixed"]) return false; else return true;
  }
  else
  // smax reached
  if(s>=P_int["smax"]) return false;
  else
  // std
  if(R<P_double["eps_s"] && s>=P_int["smin"]) return false; else return true;
}


signal_analyzer::signal_analyzer(TExperiment* _ex, double _relax, double _min_diff_threshold) : TExperiment_ref(_ex), relax(_relax), min_diff_threshold(_min_diff_threshold)
{
  count=0;
  t=0.0;

  uc=0.0;
  descending=false;
  min_t=0.0;
  threshold_reached=false;
  min_interval=0.0;
}

void signal_analyzer::step(double _t, double u)
{
  t_prev=t;
  uc_prev=uc;

  t=_t;
  uc=relax*u+(1.0-relax)*uc_prev;
  if(uc<=uc_prev)
  {
    descending=true;
  }
  else if(descending)
  {
    min_handler();
    descending=false;
  }
}

void signal_analyzer::min_handler()
{
  count++;
  // update values
  double min_t_old=min_t;
  double min_value_old=min_value;
  min_t=t_prev;
  min_interval=min_t-min_t_old;
  min_value=uc_prev;
  min_value_diff_rel=abs(min_value-min_value_old)/abs(0.5*(min_value_old+min_value));

  flog<<"Signal_analyzer: minimum found t="<<min_t<<" value="<<min_value<<" interval="<<min_interval<<" diff="<<min_value_diff_rel<<endl;


  if(min_value_diff_rel<min_diff_threshold && count>10)
  {
    min_diff_threshold_handler();
  }
}
void signal_analyzer::min_diff_threshold_handler()
{
  threshold_reached=true;
  flog<<"Signal_analyzer: threshold reached"<<endl;
}
