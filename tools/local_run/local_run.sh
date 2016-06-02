#!/bin/sh

if [[ -f hydro3.lock ]] ; then 
  echo "Error hydro3.lock exists. No action."
  exit
fi
./hydro start.hydroconf 1> local_output.txt 2> local_error.txt & echo $! > hydro3.pid

