/*******************************************************************
*****************    Parameters Variator (PVar)     ****************
********************************************************************/

#include "pvar.hpp"

PVar::PVar(TConsole* parent)
{
  console=parent;
  running=false;
}

void PVar::new_variator(string arg)
{
  VL.insert_to_end(PVar_single(arg));
}

void PVar::init()
{
  VL.change_N(0);
}

void PVar::start()
{
  thread_ptr = std::make_shared<std::thread>(&PVar::thread, this);
  running=true;
}

void PVar::thread()
{
  variations_count=0;
  console->CP_int.set("pvar_count", variations_count);
  cout<<endl<<"PVar started"<<endl;
  step(0);
  cout<<endl<<"PVar finished"<<endl;
  running=false;
}

void PVar::step(int depth)
{
  if(depth<VL.N)
  {
    PVar_single& V=VL[depth];
    for(int k=0; k<V.values.N; k++)
    {
      V.current=V.values[k];
      step(depth+1);
    }
  }
  else
  {
    variations_count++;
    console->CP_int["pvar_count"]=variations_count;
    console->cflog<<"pvar "<<variations_count<<endl;

    // wait for termination of current experiments
    while(console->pending_count>0)
    {
      sleep(1);
    }

    if(console->CP_string("pvar_script"))
    {
      try
      {
        console->cmd_run(console->CP_string["pvar_script"]);
      }
      catch(string msg)
      {
        console->logger_error() << "pvar: " + msg;
      }
    }
    else
    {
      cout<<"Warning: pvar_script is not defined"<<endl;
    }
  }
}

void PVar::term()
{
  thread_ptr->join();
  table_close();
}

void PVar::table_open(string arg)
{
  stringstream buf;
  buf.str(arg);

  string filename;
  buf>>filename;
  console->open_file(ftable,filename);
  console->cflog<<"Opened pvar_table file: '"<<filename<<"'"<<endl;

  table_columns.change_N(0);
  while(!buf.eof())
  {
    string pname;
    buf>>pname;
    table_columns.insert_to_end(pname);
  }

  for(int i=0; i<table_columns.N; i++)
  {
    ftable<<table_columns[i]<<";";
  }
  ftable<<endl;
}

void PVar::table_record(TExperiment* ex)
{
  for(int i=0; i<table_columns.N; i++)
  {
    ftable<<console->get_parameter_value(table_columns[i],ex)<<";";
  }
  ftable<<endl;
}

void PVar::table_close()
{
  ftable.close();
}


PVar_single::PVar_single()
{
  values.change_N(0);
}

PVar_single::PVar_single(string arg)
{
  stringstream buf;
  buf.str(arg);
  buf>>ptype;
  buf>>pname;
  buf>>vtype;
  if(vtype=="step_add")
  {
    nv_step_add(buf);
  }
  else
  if(vtype=="step_mul")
  {
    nv_step_mul(buf);
  }
  else
  if(vtype=="range_add")
  {
    nv_range_add(buf);
  }
  else
  if(vtype=="range_mul")
  {
    nv_range_mul(buf);
  }
  else
  if(vtype=="enum" || vtype=="enumeration")
  {
    nv_enum(buf);
  }
  else
  {
    throw string("Unknown variator type: '"+vtype+"'");
  }
}

void PVar_single::nv_step_add(stringstream& buf)
{
  if(ptype=="int")
  {
    double p0, eps;
    int N;
    buf>>p0>>eps>>N;
    values.change_N(N);
    double p=p0;
    for(int k=0; k<N; k++)
    {
      values[k]=convert_to_string(myround(p));
      p+=eps;
    }
  }
  else
  if(ptype=="double")
  {
    double p0, eps;
    int N;
    buf>>p0>>eps>>N;
    values.change_N(N);
    double p=p0;
    for(int k=0; k<N; k++)
    {
      values[k]=convert_to_string(p);
      p+=eps;
    }
  }
  else
  if(ptype=="vect")
  {
    column<double> p0, eps;
    int N;
    buf>>p0>>eps>>N;
    values.change_N(N);
    column<double> p=p0;
    for(int k=0; k<N; k++)
    {
      values[k]=convert_to_string(p);
      p=p+eps;
    }
  }
  else
  {
    throw string("ptype="+ptype+" is incompatible with vtype="+vtype);
  }
}

void PVar_single::nv_step_mul(stringstream& buf)
{
  if(ptype=="int")
  {
    double p0, eps;
    int N;
    buf>>p0>>eps>>N;
    values.change_N(N);
    double p=p0;
    for(int k=0; k<N; k++)
    {
      values[k]=convert_to_string(myround(p));
      p*=eps;
    }
  }
  else
  if(ptype=="double")
  {
    double p0, eps;
    int N;
    buf>>p0>>eps>>N;
    values.change_N(N);
    double p=p0;
    for(int k=0; k<N; k++)
    {
      values[k]=convert_to_string(p);
      p*=eps;
    }
  }
  else
  {
    throw string("ptype="+ptype+" is incompatible with vtype="+vtype);
  }
}

void PVar_single::nv_range_add(stringstream& buf)
{
  if(ptype=="int")
  {
    double pL, pR;
    int N;
    buf>>pL>>pR>>N;
    values.change_N(N);
    for(int k=0; k<N; k++)
    {
      double d=double(k)/(N-1);
      double p=pL*(1.0-d)+pR*d;
      values[k]=convert_to_string(myround(p));
    }
  }
  else
  if(ptype=="double")
  {
    double pL, pR;
    int N;
    buf>>pL>>pR>>N;
    values.change_N(N);
    for(int k=0; k<N; k++)
    {
      double d=double(k)/(N-1);
      double p=pL*(1.0-d)+pR*d;
      values[k]=convert_to_string(p);
    }
  }
  else
  if(ptype=="vect")
  {
    column<double> pL, pR;
    int N;
    buf>>pL>>pR>>N;
    values.change_N(N);
    for(int k=0; k<N; k++)
    {
      double d=double(k)/(N-1);
      column<double> p=pL*(1.0-d)+pR*d;
      values[k]=convert_to_string(p);
    }
  }
  else
  {
    throw string("ptype="+ptype+" is incompatible with vtype="+vtype);
  }
}

void PVar_single::nv_range_mul(stringstream& buf)
{
  if(ptype=="int")
  {
    double pL, pR;
    int N;
    buf>>pL>>pR>>N;
    values.change_N(N);
    for(int k=0; k<N; k++)
    {
      double d=double(k)/(N-1);
      double p=pow(pL,1.0-d)*pow(pR,d);
      values[k]=convert_to_string(myround(p));
    }
  }
  else
  if(ptype=="double")
  {
    double pL, pR;
    int N;
    buf>>pL>>pR>>N;
    values.change_N(N);
    for(int k=0; k<N; k++)
    {
      double d=double(k)/(N-1);
      double p=pow(pL,1.0-d)*pow(pR,d);
      values[k]=convert_to_string(p);
    }
  }
  else
  {
    throw string("ptype="+ptype+" is incompatible with vtype="+vtype);
  }
}

void PVar_single::nv_enum(stringstream& buf)
{
  int k=0;
  while(!buf.eof())
  {
    values.change_N(sup_power2(k));
    buf>>values[k];
    k++;
    if(k>1E5) throw string("pvar_enum: k>1E5");
  }
  values.change_N(k);
}
