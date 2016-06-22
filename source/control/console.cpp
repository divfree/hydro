#include "console.hpp"
#include <chrono>
#include <thread>

string CMap_name[Map_number]={"c_bool", "c_double", "c_int", "c_string"};

void sleep(int seconds) {
  std::this_thread::sleep_for(std::chrono::seconds(seconds));
}

void sleep_ms(int ms) {
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

void TConsole::cmd_exit(string) {
  scheduler_term();
  int expn=undone_exp_count();
  if(expn>0)
  {
    confirmation_yes_no("There are "+IntToStr(expn)+" undone experiments. Do you really want to exit?");
  }
  while(Experiments.length>0)
  {
    TExperiment* ex=Experiments.get_data(0);
    exp_del(ex);
  }
  flag_exit=true;
}

void TConsole::cmd_echo(string arg) {
  std::cout << arg;
}

void TConsole::cmd_add_experiment(string arg) {
  stringstream buf;
  buf.str(arg);
  string exp_name;
  buf>>exp_name;

  if(!Experiments(exp_name))
  {
    TExperiment* new_exp=new TExperiment(this);
    new_exp->name=exp_name;
    Experiments.set(exp_name, new_exp);
    cur_exp_name=exp_name;
  } else
  {
    throw string("Experiment '"+exp_name+" already exists");
  }
}

void TConsole::exp_del(TExperiment* ex) {
  cflog<<"("<<ex->name<<", "<<ex->P_string["name"]<<" ) deleted"<<endl;
  if(!ex->st_term)
  {
    exp_term(ex);
  }
  Experiments.del(ex->name);
  delete ex;
}

void TConsole::cmd_del_experiment(string arg) {
  stringstream buf;
  buf.str(arg);
  string exp_name;
  buf>>exp_name;

  if(Experiments.exist(exp_name))
  {
    TExperiment* ex=Experiments[exp_name];
    if(ecast(ex->P_bool("started")) && !ecast(ex->P_bool("thread_completed")))
    confirmation_yes_no("Do you really want to delete the experiment?");
    exp_del(ex);
    if(cur_exp_name==exp_name) cur_exp_name="";
  } else
  {
    throw string("Experiment '"+exp_name+" doesn't exist");
  }
}

void TConsole::cmd_current_experiment(string arg) {
  stringstream buf;
  buf.str(arg);
  string exp_name;
  buf>>exp_name;

  if(exp_name=="") cmd_echo(cur_exp_name);
    else if(Experiments(exp_name)==0) throw string("Experiment '"+exp_name+"' doesn't exist");
      else cur_exp_name=exp_name;
}

void TConsole::cmd_init(string) {
  TExperiment &cur_exp=*check_cur_exp();
  cur_exp.init();
}

void TConsole::cmd_status(string) {
  if(cur_exp_name!="")
  {
    TExperiment &cur_exp=*check_cur_exp();
    cout<<"Experiment '"<<cur_exp_name<<"' status: "<<cur_exp.status_name()<<endl;
  } else
  {
    cout<<"No current experiment"<<endl;
  }
  cout<<"Total experiments number: "<<Experiments.length<<endl;
  cout<<"Virtual memory usage: "<<sysinfo::virtual_usage_kb()/1000<<"MB"<<endl;
  cout<<"Physical memory usage: "<<sysinfo::physical_usage_kb()/1000<<"MB"<<endl;
  cout<<"Threads: "<<sysinfo::threads()<<endl;
}

void TConsole::cmd_wait_for_experiments_completion(string) {
  while(undone_exp_count()>0 || pvar->running)
  {
    sleep(1);
  }
}

void TConsole::cmd_run(string arg) {
  stringstream buf;
  buf.str(arg);
  string filename="";
  filename=get_string_with_ws(buf);
  ifstream fin(filename);

  if(!fin.good()) throw string("File '"+filename+"' is not good");
  while(fin.good())
  {
    string str_input="";
    getline(fin, str_input);
    str_input=pname_parser(str_input);
    string cur_name, cur_arg;
    extract_cmd_name_and_arg(str_input, cur_name, cur_arg);
    TCmd* ptr=Commands(cur_name);
    if(ptr)
    {
      (this->**ptr)(cur_arg);
    }
    else throw string("Unknown command: '"+cur_name+"'");
  }
  fin.close();
}

void TConsole::exp_start(TExperiment* experiment) {
  TExperiment &cur_exp=*experiment;
  if(cur_exp.st_thread)
  {
    cout<<"Computation in experiment '"+cur_exp_name+"' is already started"<<endl;
  } else
  if(cur_exp.st_init)
  {
    cur_exp.st_pending=true;
    pending_count++;
    cur_exp.P_bool.set("started", true);
  } else
  {
    throw string("Experiment '"+cur_exp_name+"' is not initialized");
  }
}

void TConsole::exp_create_thread(TExperiment* experiment) {
  TExperiment &cur_exp=*experiment;
  if(cur_exp.st_pending)
  {
    cur_exp.thread_ptr =
        std::make_shared<std::thread>(&TExperiment::thread, &cur_exp);
    cur_exp.st_pending=false;
    cur_exp.st_thread=true;
    pending_count--;
    threads_count++;
    cout<<"("<<cur_exp.name<<") started"<<endl;
    cflog<<"("<<cur_exp.name<<") started"<<endl;
  }
  else
  {
    throw string("Experiment '"+cur_exp_name+"' is not pending");
  }
}
void TConsole::cmd_start(string) {
  exp_start(check_cur_exp());
}

void TConsole::cmd_start_in_main_thread(string) {
  TExperiment &cur_exp=*check_cur_exp();
  if(ecast(cur_exp.P_bool("started")))
  {
    cout<<"Computation in experiment '"+cur_exp_name+"' is already started"<<endl;
  } else
  if(cur_exp.st_init)
  {
    cout<<"("<<cur_exp.name<<") started"<<endl;
    cur_exp.st_thread=true;
    cur_exp.thread();
  } else
  {
    throw string("Experiment '"+cur_exp_name+"' is not initialized");
  }
}

void TConsole::exp_stop(TExperiment* ex) {
  if(ex->st_thread)
  {
    ex->P_bool["terminate"]=true;
    ex->sig_term=true;

    cout<<"Waiting for thread termination...";
    ex->thread_ptr->join();
    ex->thread_ptr.reset();
    cout<<"OK"<<endl;
  } else
  {
    cout<<"Computation in experiment '"+ex->name+"' is not started"<<endl;
  }
}

void TConsole::cmd_stop(string) {
  exp_stop(check_cur_exp());
}

void TConsole::exp_term(TExperiment* ex) {
  if(ex->st_init)
  {
    exp_stop(ex);
    ex->term();
    cout<<"Experiment '"<<ex->name<<"' successfully terminated"<<endl;
  } else
  {
    cout<<"Experiment '"<<ex->name<<"' is not initialized"<<endl;
  }
}

void TConsole::cmd_term(string) {
  exp_term(check_cur_exp());
}

void TConsole::cmd_list_exp(string) {
  cout<<"Experiments: "<<Experiments.length<<endl;
  for(int i=0; i<Experiments.length; i++) cout<<Experiments.get_key(i)<<" ";
  cout<<endl;
  cout<<endl;

  int running=0;
  for(int i=0; i<Experiments.length; i++)
    if(Experiments.get_data(i)->st_thread) running++;
  cout<<"Running experiments: "<<running<<endl;
  if(running)
  {
    for(int i=0; i<Experiments.length; i++)
      if(Experiments.get_data(i)->st_thread) cout<<Experiments.get_key(i)<<" ";
    cout<<endl;
  }
  cout<<endl;

  int done=0;
  for(int i=0; i<Experiments.length; i++)
    if(Experiments.get_data(i)->st_finished) done++;
  cout<<"Done experiments: "<<done<<endl;
  if(done)
  {
    for(int i=0; i<Experiments.length; i++)
      if(Experiments.get_data(i)->st_finished) cout<<Experiments.get_key(i)<<" ";
    cout<<endl;
  }
  cout<<endl;
}

TExperiment* TConsole::check_cur_exp() {
  if(!Experiments.exist(cur_exp_name))
  {
    throw string("Experiment '"+cur_exp_name+"' doesn't exist");
  }
  return Experiments[cur_exp_name];
}

// TODO: Parameters parsed on evaluation
void TConsole::cmd_set(string arg) {
  stringstream buf;
  buf.str(arg);
  string p_type;
  buf>>p_type;

  bool found=false;

  for(int k=0; k<CMap_number; k++)
  if(CMap_name[k]==p_type)
  {
    found=true;
    string name;
    buf>>name;
    bool running_smth=false;
    for(int i=0; i<Experiments.length; i++)
    {
      running_smth=running_smth || Experiments.get_data(i)->st_thread;
    }
    if(running_smth && !CMap[k]->exist(name)) throw string("Can't add new console parameter during computation of any experiment");
    CMap[k]->read_data_by_key(buf, name);
  }

  if(found) return;

  TExperiment& ex=*check_cur_exp();

  for(int k=0; k<Map_number; k++)
  if(Map_name[k]==p_type)
  {
    found=true;
    string name;
    buf>>name;
    if(ex.st_thread && !ex.Map[k]->exist(name)) throw string("Can't add new parameter during computation");
    //if(!ex.Map[k]->PK.exist(name)) cout<<"Warning: '"<<Map_name[k]<<" "<<name<<"' is not known parameter"<<endl;
    ex.Map[k]->read_data_by_key(buf, name);
  }
  if(!found) throw string("Unknown type '"+p_type+"'");
}

void TConsole::cmd_del(string arg) {
  stringstream buf;
  buf.str(arg);
  string name;
  buf>>name;
  bool found=false;
  for(int k=0; k<CMap_number; k++)
  {
    if(CMap[k]->del(name))
    {
      cout<<"Deleted "<<Map_name[k]<<" "<<name<<endl;
      found=true;
    }
  }
  if(found) return;

  TExperiment& ex=*check_cur_exp();
  for(int k=0; k<Map_number; k++)
  {
    if(ex.Map[k]->del(name))
    {
      cout<<"Deleted "<<Map_name[k]<<" "<<name<<endl;
      found=true;
    }
  }
  if(!found) cout<<"Parameter '"+name+"' not found";
}

void TConsole::cmd_list_parameters(string arg) {
  stringstream buf;
  buf.str(arg);
  string type;
  string s1="", s2="";
  buf>>s1>>s2;

  bool names_only=false;
  if(s1=="-n")
  {
    names_only=true;
    type=s2;
  } else
  if(s2=="-n")
  {
    names_only=true;
    type=s1;
  } else type=s1;

  if(type=="" && names_only)
  {
    cout<<"Console"<<endl;
    for(int k=0; k<CMap_number; k++)
    {
      if(CMap[k]->get_length()>0)
      {
        if(k) cout<<endl;
        cout<<"TYPE "<<CMap_name[k]<<endl;
        for(int i=0; i<CMap[k]->get_length(); i++)
        {
          cout<<(i?", ":"")<<CMap[k]->get_key(i);
        }
        cout<<endl;
      }
    }

    cout<<endl;

    if(Experiments(cur_exp_name))
    {
      TExperiment& ex=*check_cur_exp();
      cout<<"Experiment: "<<cur_exp_name<<endl;
      for(int k=0; k<Map_number; k++)
      {
        if(ex.Map[k]->get_length()>0)
        {
          if(k) cout<<endl;
          cout<<"TYPE "<<Map_name[k]<<endl;
          for(int i=0; i<ex.Map[k]->get_length(); i++)
          {
            cout<<(i?", ":"")<<ex.Map[k]->get_key(i);
          }
          cout<<endl;
        }
      }
    }
    return;
  }


  for(int k=0; k<CMap_number; k++)
  {
    if(CMap_name[k]==type)
    {
      cout<<"TYPE "<<CMap_name[k]<<" (Console)"<<endl;
      for(int i=0; i<CMap[k]->get_length(); i++)
      {
        if(names_only)
        {
          cout<<(i?", ":"")<<CMap[k]->get_key(i);
        } else
        {
          cout<<CMap[k]->get_key(i)<<"=";
          CMap[k]->write_data_by_index(cout, i);
          cout<<endl;
        }
      }
      cout<<endl;
    }

    if(CMap[k]->exist(type))
    {
      cout<<CMap_name[k]<<" "<<type<<"=";
      CMap[k]->write_data_by_key(cout, type);
      cout<<endl;
    }
  }

  if(Experiments(cur_exp_name))
  {
    TExperiment& ex=*check_cur_exp();
    for(int k=0; k<Map_number; k++)
    {
      if(Map_name[k]==type)
      {
        cout<<"TYPE "<<Map_name[k]<<" (Exp)"<<endl;
        for(int i=0; i<ex.Map[k]->get_length(); i++)
        {
          if(names_only)
          {
            cout<<(i?", ":"")<<ex.Map[k]->get_key(i);
          } else
          {
            cout<<ex.Map[k]->get_key(i)<<"=";
            ex.Map[k]->write_data_by_index(cout, i);
            cout<<endl;
          }
        }
        cout<<endl;
      }
      if(ex.Map[k]->exist(type))
      {
        cout<<Map_name[k]<<" "<<type<<"=";
        ex.Map[k]->write_data_by_key(cout, type);
        cout<<endl;
      }
    }
  }
}

string TConsole::get_parameter_value(string arg, TExperiment* ex) {
  stringstream buf;
  buf.str(arg);
  string name;
  buf>>name;
  string res;
  bool found=false;

  stringstream buf2;
  buf.str("");

  for(int k=0; k<CMap_number; k++)
  {
    if(CMap[k]->exist(name))
    {
      if(!found)
      {
        CMap[k]->write_data_by_key(buf2, name);
        buf2>>res;
        found=true;
      }
      else
      {
        throw string("More than one parameter with name '"+name+"'");
      }
    }
  }

  if(ex || Experiments(cur_exp_name))
  {
    if(!ex) ex=check_cur_exp();
    for(int k=0; k<Map_number; k++)
    {
      if(ex->Map[k]->exist(name))
      {
        if(!found)
        {
          ex->Map[k]->write_data_by_key(buf2, name);
          buf2>>res;
          found=true;
        }
        else
        {
          throw string("More than one parameter with name '"+name+"'");
        }
      }
    }
  }
  if(found) return res; else throw string("Parameter '"+name+"' not found");
}

void TConsole::cmd_value(string arg) {
  cout<<get_parameter_value(arg);
}


void TConsole::cmd_pvar_new(string arg) {
  pvar->new_variator(arg);
}

void TConsole::cmd_pvar_start(string) {
  pvar->start();
}

void TConsole::cmd_pvar_set(string) {
  for(int i=0; i<pvar->VL.N; i++)
  {
    PVar_single& V=pvar->VL[i];
    cmd_set(V.ptype+" "+V.pname+" "+V.current);
  }
}

void TConsole::cmd_pvar_table(string arg) {
  pvar->table_open(arg);
}

void TConsole::cmd_open_log(string arg) {
  stringstream buf;
  buf.str(arg);
  string filename;
  buf>>filename;
  cflog.open(filename);
  cflog<<"Console log"<<endl;
}

void TConsole::cmd_test(string) {
}

void TConsole::cmd_test2(string) {

}

void TConsole::cmd_cd(string arg) {
  stringstream buf;
  buf.str(arg);
  string dir=get_string_with_ws(buf);
  if(dir=="")
  {
    cout<<get_current_directory()<<endl;
  } else
  {
    set_current_directory(dir);
  }
}

void TConsole::cmd_dir(string arg) {
  using namespace boost::filesystem;
  stringstream buf;
  buf.str(arg);
  string dir=get_string_with_ws(buf);
  if(dir=="") dir=".";
  try
  {
    boost::filesystem::path p(dir.c_str());
    copy(directory_iterator(p), directory_iterator(),
            std::ostream_iterator<directory_entry>(cout, "\n"));
  }
  catch (const filesystem_error& ex)
  {
    cout << ex.what() << '\n';
  }
}

void TConsole::cmd_sleep(string arg) {
  stringstream buf;
  buf.str(arg);
  int n;
  buf>>n;
  sleep(n);
}

void TConsole::cmd_usleep(string arg) {
  stringstream buf;
  buf.str(arg);
  int n;
  buf>>n;
  sleep_ms(n);
}

void TConsole::cmd_nop(string) {

}

void TConsole::cmd_hello(string) {
  cout<<"Hi! Nice to see you."<<endl;
}

void TConsole::cmd_save(string arg) {
  TExperiment* ex=check_cur_exp();
  stringstream buf;
  buf.str(arg);
  string filename="";
  filename=get_string_with_ws(buf);

  ex->module->save_state(filename);

  cout<<"Context of '"<<cur_exp_name<<"' successfully saved"<<endl;
}

void TConsole::cmd_load(string arg) {
  TExperiment* ex=check_cur_exp();

  stringstream buf;
  buf.str(arg);
  string filename="";
  filename=get_string_with_ws(buf);

  ex->module->load_state(filename);

  cout<<"Context of '"<<cur_exp_name<<"' successfully loaded"<<endl;
}

void TConsole::cmd_help(string) {
  cout<<"Commands list:"<<endl;
  for(int i=0; i<Commands.length; i++)
  {
    string name=Commands.get_key(i);
    string desc;
    string* ptr=Commands_desc(name);
    if(ptr) desc=*ptr;
    cout<<name<<" - "<<desc<<endl;
  }
}

void TConsole::cmd_mkdir(string arg) {
  stringstream buf;
  buf.str(arg);
  string dir=get_string_with_ws(buf);
  create_directory(dir);
}


TConsole::TConsole()
    : logger_error("Error: "),
      logger_warning("Warning: "),
      logger_info(),
      logger_prompt(true) {
  finished=false;

  TParameter<string>* CMap1[CMap_number]={&CP_bool, &CP_double, &CP_int, &CP_string};
  for(int k=0; k<CMap_number; k++)
  {
    CMap[k]=CMap1[k];
  }

  // TODO: Automatic console command registration

  flag_exit=false;
  Commands.set("exit",&TConsole::cmd_exit);
  Commands.set("echo",&TConsole::cmd_echo);
  Commands.set("add_exp", &TConsole::cmd_add_experiment);
  Commands.set("ae", &TConsole::cmd_add_experiment);
  Commands.set("del_exp", &TConsole::cmd_del_experiment);
  Commands.set("de", &TConsole::cmd_del_experiment);
  Commands.set("init", &TConsole::cmd_init);
  Commands.set("status", &TConsole::cmd_status);
  Commands.set("start", &TConsole::cmd_start);
  Commands.set("start_in_main_thread", &TConsole::cmd_start_in_main_thread);
  Commands.set("stop", &TConsole::cmd_stop);
  Commands.set("cur_exp",&TConsole::cmd_current_experiment);
  Commands.set("ce",&TConsole::cmd_current_experiment);
  Commands.set("list_exp",&TConsole::cmd_list_exp);
  Commands.set("le",&TConsole::cmd_list_exp);
  Commands.set("set",&TConsole::cmd_set);
  Commands.set("list_parameters",&TConsole::cmd_list_parameters);
  Commands.set("lp",&TConsole::cmd_list_parameters);
  Commands.set("value",&TConsole::cmd_value);
  Commands.set("del",&TConsole::cmd_del);
  Commands.set("help",&TConsole::cmd_help);
  Commands.set("run",&TConsole::cmd_run);
  Commands.set("sleep",&TConsole::cmd_sleep);
  Commands.set("usleep",&TConsole::cmd_usleep);
  Commands.set("nop",&TConsole::cmd_nop);
  Commands.set("#",&TConsole::cmd_nop);
  Commands.set("",&TConsole::cmd_nop);
  Commands.set("cd",&TConsole::cmd_cd);
  Commands.set("dir",&TConsole::cmd_dir);
  Commands.set("save",&TConsole::cmd_save);
  Commands.set("load",&TConsole::cmd_load);
  Commands.set("term",&TConsole::cmd_term);
  Commands.set("test",&TConsole::cmd_test);
  Commands.set("hello",&TConsole::cmd_hello);
  Commands.set("mkdir",&TConsole::cmd_mkdir);
  Commands.set("wait_for_completion",&TConsole::cmd_wait_for_experiments_completion);
  Commands.set("pvar_new",&TConsole::cmd_pvar_new);
  Commands.set("pvar_start",&TConsole::cmd_pvar_start);
  Commands.set("pvar_set",&TConsole::cmd_pvar_set);
  Commands.set("pvar_table",&TConsole::cmd_pvar_table);
  Commands.set("open_log",&TConsole::cmd_open_log);

  Commands_desc.set("help", "view this help");
  Commands_desc.set("lp", "alias for list_parameters");
  Commands_desc.set("le","alias for list_exp");
  Commands_desc.set("set","create or modify parameter");
  Commands_desc.set("del","delete parameter");
  Commands_desc.set("cur_exp","change current experiment");
  Commands_desc.set("ce","alias for cur_exp");
  Commands_desc.set("init","initialize current experiment using known parameters");
  Commands_desc.set("exit","terminate program");
  Commands_desc.set("start","start the computation in current experiment");
  Commands_desc.set("stop","stop the computation in current experiment");
  Commands_desc.set("add_exp","add new experiment");
  Commands_desc.set("ae","alias for add_exp");
  Commands_desc.set("term","terminate the experiment");
  Commands_desc.set("run","run the script");
  Commands_desc.set("sleep","wait for <n> seconds");
  Commands_desc.set("usleep","wait for <n> milliseconds");
  Commands_desc.set("copy_comp","alias for copy_U_comp_to_grid_double");
  Commands_desc.set("copy_U_comp_to_grid_double","syntax example: u grid_u");
  Commands_desc.set("wait_for_completion","wait for all experiments termination");
  Commands_desc.set("value","provides the value of parameter");
  Commands_desc.set("pvar_table","<filename> <p1> <p2> ...");

  pvar=new PVar(this);

  scheduler_init();
}

TConsole::~TConsole() {
  if(!finished)
  {
    cmd_exit("");
  }
  cflog<<"Console termination"<<endl;
  cflog.close();
  delete pvar;
}

void extract_cmd_name_and_arg(string str, string& name, string& arg) {
  stringstream buf;
  buf.str(str);
  buf>>name;
  if(name.length()>0 && name[0]=='#')
  {
    name="#";
    arg="";
  } else
  {
    getline(buf, arg);
  }
}

string TConsole::prompt() {
  //string str=get_time("%X");
  string str=cur_exp_name;
  return get_current_directory()+":"+str+"> ";
}

void TConsole::execute() {
  while (!flag_exit) {
    try {
      cout<<prompt();
      string str_input="";
      getline(cin, str_input);
      str_input=pname_parser(str_input);
      string name, arg;
      extract_cmd_name_and_arg(str_input, name, arg);
      TCmd* ptr=Commands(name);
      if(ptr) {
        (this->**ptr)(arg);
      } else {
        logger_error() << "Unknown command: '" << name << "'";
      }
      cout<<endl;
    }
    catch(string msg) {
      logger_error() << msg;
    }
  }
  finished=true;
}

void TConsole::check_directory(string dirname) {
  if(!directory_exists(dirname))
  {
    char c=0;
    while(c!='y' && c!='n')
    {
      logger_prompt() << "Directory '"+dirname+"' doesn't exist. Create? [y/n] ";
      char c;
      if(ecast(CP_bool("force_yes")))
      {
        c='y';
      } else
      {
        cin>>c;
      }
      if(c=='n')
      {
        break;
      } else
      if(c=='y')
      {
        create_directory(dirname);
        break;
      }
    }
  }
}

void TConsole::check_file_good(ofstream& file, string filename) {
  if(!file.good())
  {
    throw string("Can't open file '"+filename+"'");
  }
}

void TConsole::check_file_good(ifstream& file, string filename) {
  if(!file.good())
  {
    throw string("Can't open file '"+filename+"'");
  }
}

void TConsole::confirmation_yes_no(string msg) {
  while(true)
  {
    logger_prompt() << msg << " [y/n] ";
    char c;
    if (ecast(CP_bool("force_yes"))) {
      c='y';
    } else {
      cin>>c;
    }
    if (c == 'n') {
      throw string("Operation canceled");
    } else if (c == 'y') {
      break;
    }
  }
}

void TConsole::check_file_exists(string filename) {
  if(file_exists(filename) && ! ecast(CP_bool("force_overwrite")))
  {
    while(true)
    {
      logger_prompt() << "File '" << filename
          << "' already exists. Overwrite? [y/n/a(all)] ";
      char c;
      if(ecast(CP_bool("force_yes")))
      {
        c='y';
      } else
      {
        cin>>c;
      }
      if(c=='n')
      {
        throw string("Overwriting of existing file canceled");
      } else
      if(c=='y') break; else
      if(c=='a')
      {
        CP_bool.set("force_overwrite", true);
        break;
      }
    }
  }
}

void TConsole::open_file(ofstream& fout, string filename) {
  string path=filename;
  check_file_exists(path);
  fout.open(path);
  check_file_good(fout,filename);
  logger_info() << "File '" << filename << "' is opened (write-mode)";
  fout.precision(16);
}

void TConsole::open_file(ifstream& fin, string filename) {
  string path = filename;
  fin.open(path);
  check_file_good(fin, filename);
  logger_info() << "File '" << filename << "' is opened (read-mode)";
}

int TConsole::running_exp_count() {
  int count=0;
  for(int i=0; i<Experiments.length; i++)
  {
    if(Experiments.get_data(i)->st_thread) count++;
  }
  return count;
}


int TConsole::undone_exp_count() {
  int count=0;
  for(int i=0; i<Experiments.length; i++)
  {
    TExperiment* ex=Experiments.get_data(i);
    if(ex->st_init && !ex->st_term && !ex->st_error)  count++;
  }
  return count;
}

int TConsole::available_threads_number() {
  if(CP_int("max_threads_number")) return CP_int["max_threads_number"]; else return 4;
}

void TConsole::scheduler_init() {
  threads_count=0;
  pending_count=0;
  scheduler_terminate=false;
  scheduler_terminate_done=false;
  scheduler_thread_ptr =
      std::make_shared<std::thread>(&TConsole::scheduler_thread, this);
}

void TConsole::scheduler_thread() {
  while(!scheduler_terminate)
  {
    /*for(int i=0; i<Experiments.length; i++)
    {
      TExperiment* ex=Experiments.get_data(i);
      if(ex->st_init && !ex->st_thread && (ex->st_term || ex->st_error))
      {
        cflog<<"Removing of "<<ex->name<<endl;
        exp_del(Experiments.get_data(i));
        i--;
      }
    }*/

    if(threads_count<available_threads_number() && pending_count>0)
    {
      int k=0;
      while(threads_count<available_threads_number() && pending_count>0 && k<Experiments.length)
      {
        TExperiment* ex=Experiments.get_data(k);
        if(ex->st_pending)
        {
          exp_create_thread(ex);
        }
        k++;
      }
    }
    sleep_ms(100);
  }
  scheduler_terminate_done=true;
}

void TConsole::scheduler_term() {
  cout<<endl<<"Scheduler termination...";
  scheduler_terminate=true;
  scheduler_thread_ptr->join();
  cout<<"done"<<endl;
}

bool TConsole::valid_char_for_parameter_name(char c) {
  return (c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9') || (c=='_');
}

/*********************************************
   $(pname) and $pname constructions parser
INPUT: command line
OUTPUT: command line with all $(pname) and $pname
replaced by the value of parameter.
The value is similar to cmd_value output
********************************************/
string TConsole::pname_parser(string cmdline) {
  size_t L=cmdline.length();
  if(L==0 || cmdline[0]=='#') return cmdline;
  string res="";
  const string& s=cmdline;
  int state=0;
  size_t k=0;
  string P;
  while(state>=0)
  {
    switch(state)
    {
    case 0:
      if(k>=L) { state=-1; break; }
      if(s[k]=='$') state=2; else state=1;
      break;
    case 1:
      res+=s[k]; k++;
      if(k>=L) { state=-1; break; }
      state=0;
      break;
    case 2:
      k++; P="";
      if(k>=L) { state=-2; break; }
      if(s[k]=='(') state=3;
      else if(valid_char_for_parameter_name(s[k])) state=7;
      else state=-4;
      break;
    case 3:
      k++;
      if(k>=L) state=-3;
      else if(!valid_char_for_parameter_name(s[k])) state=-4;
      else state=4;
      break;
    case 4:
      P+=s[k];
      k++;
      if(k>=L) state=-2;
      else if(valid_char_for_parameter_name(s[k])) state=4;
      else if(s[k]==')') state=5;
      else state=-4;
      break;
    case 5:
      k++;
      state=8;
      break;
    case 7:
      P+=s[k];
      k++;
      if(k>=L || !valid_char_for_parameter_name(s[k])) state=8;
      else state=7;
      break;
    case 8:
      res+=get_parameter_value(P);
      state=0;
      break;
    default:
      throw string("parser: Unknown state: "+IntToStr(state));
    }
  }
  string last_char(1,s[k]);
  if(state==-2) throw string("parser: Invalid line end while expecting a parameter name");
  else if(state==-3) throw string("parser: Matching parentheses ')' not found");
  else if(state==-4) throw string("parser: Invalid char '"+last_char+"' in parameter name");
  return res;
}
