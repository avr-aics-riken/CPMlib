#! /bin/sh
#
# at .bashrc
#
# Compiler options:
#
#   --with-comp=INTEL|FJ;      If compiler does not fall under the category, this option is blank.
#   --with-ompi=/hoge;         In case of using wrapper compiler, this option may be omitted.
#   --with-f90real=4|8;        Specify real type in fortran
#   --with-f90example=yes|no;  Specify compilation of fortran sample included. 
#   --host=hostname;           Specify in case of cross-compilation.
#
#
./configure --prefix=$1 \
            --with-comp=INTEL \
            --with-ompi=/opt/openmpi \
	    --with-parser=/usr/local/textparser \
            --with-f90example=no \
            CXX=icpc \
            CXXFLAGS=-O3 \
	    FC=ifort
