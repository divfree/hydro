#include "control/console.hpp"
#include <iostream>
#include <string>
#include "revision.hpp"

#ifdef MPI_ENABLE
#include <mpi.h>
#endif

int main(int argc, char* argv[])
{
  std::string lockfilename="hydro3.lock";
#ifndef MPI_ENABLE
  if (file_exists(lockfilename)) {
    std::cout << "Error: " + lockfilename +
        " exists. Another instance seems to be running." << std::endl;
    return -1;
  }
#endif

#ifdef MPI_ENABLE
  MPI_Init(NULL, NULL);
#endif

  std::ofstream lockfile(lockfilename);
  lockfile.close();

  int exitcode=-1;

  #ifdef _DEBUG
    cout<<"Warning: _DEBUG defined"<<endl;
  #endif
  #ifdef _DEBUG2
    cout<<"Warning: _DEBUG2 defined"<<endl;
  #endif
  #ifdef _DEBUG3
    cout<<"Warning: _DEBUG3 defined"<<endl;
  #endif

  TConsole console;

  if(argc>2)
  {
    cout<<"Too much arguments: "<<argc<<endl;
    exitcode=-1;
  }
  else if(argc==2) // start the specified script
  {
    try
    {
      cout<<"Starting the script..."<<endl;
      // run the script
      console.cmd_run(argv[1]);

      cout<<"Waiting for all experiments termination (may take much time)..."<<endl;
      // wait for all experiments termination
      console.cmd_wait_for_experiments_completion("");
    }
    catch(string msg)
    {
      console.error(msg);
      exitcode=-1;
    }
    cout<<endl;
    exitcode=0;
  }
  else if(argc==1 || argc==0) // interactive mode
  {
    std::cout << "Welcome to hydro3\n";
    std::cout << "revision " << kGitRevision << "\n";
    std::cout << "Type in a command" << std::endl;
    console.execute();
    cout<<endl;
    exitcode=0;
  }
  else
  {
    cout<<"Something wrong with parameters count: "<<argc<<endl;
    exitcode=-1;
  }

  if (file_exists(lockfilename)) {
    remove_file(lockfilename);
  }
  
//#ifdef MPI_ENABLE
//  MPI_Init(NULL, NULL);
//
//  int world_size;
//  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//  int world_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//
//  char processor_name[MPI_MAX_PROCESSOR_NAME];
//  int name_len;
//  MPI_Get_processor_name(processor_name, &name_len);
//
//  printf("Hello world from processor %s, rank %d"
//         " out of %d processors\n",
//         processor_name, world_rank, world_size);
//#endif

#ifdef MPI_ENABLE
  MPI_Finalize();
#endif
  
  return exitcode;
}
