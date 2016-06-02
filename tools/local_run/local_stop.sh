#!/bin/sh
if [[ ! -f hydro3.lock || -z `pidof hydro3` ]] ; then
  echo "Error: hydro3 not running. No action."
  exit
fi

if [[ ! -f hydro3.pid ]] ; then
  echo "Error: hydro3.pid not found. No action."
  exit
fi

kill `cat hydro3.pid`
rm 'hydro3.pid' 'hydro3.lock'