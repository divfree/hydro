#include "control/console.hpp"
#include <iostream>
#include <string>
#include "revision.hpp"
#include "control/logger.hpp"

#ifdef MPI_ENABLE
#include <mpi.h>
#endif

int main(int argc, char* argv[])
{
  logger::Logger logger_error("Error: ");
  logger::Logger logger_warning("Warning: ");
  logger::Logger logger_info;

  std::string lockfilename="hydro3.lock";
#ifndef MPI_ENABLE
  if (file_exists(lockfilename)) {
    logger_error() << lockfilename
        << " exists. Another instance seems to be running.";
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
  logger_warning() << "_DEBUG defined";
#endif
#ifdef _DEBUG2
  logger_warning() << "_DEBUG2 defined";
#endif
#ifdef _DEBUG3
  logger_warning() << "_DEBUG3 defined";
#endif

  TConsole console;

  if(argc==2) { // start the specified script
    try
    {
      logger_info() << "Starting the script...";
      // run the script
      console.cmd_run(argv[1]);

      logger_info() << "Waiting for all experiments termination...";
      // wait for all experiments termination
      console.cmd_wait_for_experiments_completion("");
    }
    catch(string msg) {
      console.logger_error() << msg;
      exitcode=-1;
    }
    cout<<endl;
    exitcode=0;
  } else if(argc==1 || argc==0) { // interactive mode
    std::cout << "Welcome to hydro3\n";
    std::cout << "revision " << kGitRevision << "\n";
    std::cout << "Type in a command" << std::endl;
    console.execute();
    cout<<endl;
    exitcode=0;
  } else {
    logger_info() << "Incorrect number of arguments: " << argc;
    exitcode=-1;
  }

  if (file_exists(lockfilename)) {
    remove_file(lockfilename);
  }

#ifdef MPI_ENABLE
  MPI_Finalize();
#endif
  
  return exitcode;
}
