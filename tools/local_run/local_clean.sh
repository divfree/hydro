#!/bin/sh

if [[ -f hydro3.pid ]] ; then
  pid=`cat hydro3.pid`
  if ps -p $pid > /dev/null ; then
    echo "Error: hydro3 running. No action."
    exit
  fi
fi

trash=trash
mkdir $trash
rm hydro3.lock
rm hydro3.pid
rm $trash/*
mv tp1* $trash
