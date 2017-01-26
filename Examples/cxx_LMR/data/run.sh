#!/bin/sh

  rm -rf *.log

#  mpirun -np 15 ../exampleCXX_LMR leaf15.tp #leaf15.tp
  mpirun -np 4 ../exampleCXX_LMR leaf15.tp #leaf15.tp

