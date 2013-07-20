###
###################################################################################
#
# CPMlib - Cartesian Partition Manager Library
# 
# Copyright (C) 2012-2013 Institute of Industrial Science, The University of Tokyo.
# All right reserved.
#
###################################################################################
###

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([CPMlib],[1.1.2],[keno@iis.u-tokyo.ac.jp],[CPMlib])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/cpm_ParaManager.cpp])
AM_CONFIG_HEADER([config.h])

#
# revision No.
#
CPM_REVISION="20130720_1200"
AC_SUBST(CPM_REVISION)



#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/CPMlib])
if test "$prefix" = "NONE" ; then
  CPM_INST_DIR=/usr/local/CPMlib
else
  CPM_INST_DIR="$prefix"
fi



#
# platform
#
TGT_DEVICE="$($SHELL ./chk-uname)"
AC_SUBST(TGT_DEVICE)



#
# compiler
#
AC_ARG_WITH(comp, [AC_HELP_STRING([--with-comp=(INTEL|FJ)],[Specify Compiler type])], , with_comp=none)



#
# REAL environment for Fortran90
#
AC_ARG_WITH(f90real, [AC_HELP_STRING([--with-f90real=(4|8)],[Specify Fortran90 CPM_REAL size [4byte]])], , with_f90real=4)
AC_SUBST(REALOPT)
AC_SUBST(FREALOPT)

if test "$with_f90real" = "8" ; then
  REALOPT=-D_REAL_IS_DOUBLE_
  if test "$with_comp" = "INTEL" ; then
    FREALOPT="-r8"
  elif test "$with_comp" = "FJ" ; then
    FREALOPT="-CcdRR8"
  fi
else
  REALOPT=
  FREALOPT=
fi



#
# MPI setting. wrappers like mpicxx, mpicc, mpiCC, mpic++, mpiFCCpx are taking 
# into consideration, if not add in case.
#

enable_mpi="none"


CPM_CXX="$CXX"
AC_SUBST(CPM_CXX)

case "$CPM_CXX" in
  mpicxx)   enable_mpi="wrappermpi" ;;
  mpicc)    enable_mpi="wrappermpi" ;;
  mpiCC)    enable_mpi="wrappermpi" ;;
  mpic++)   enable_mpi="wrappermpi" ;;
  mpiFCCpx) enable_mpi="wrappermpi" ;;
esac



#
# MPICH environment
#
AC_SUBST(MPICH_DIR)
AC_SUBST(MPICH_CFLAGS)
AC_SUBST(MPICH_LDFLAGS)
AC_SUBST(MPICH_LIBS)
AC_SUBST(MPICH_FCLIBS)


#
# if wrapper compiler is specified, skip this part
#

if test "$enable_mpi" = "none" ; then

# MPICH
  AC_ARG_WITH(mpich, [AC_HELP_STRING([--with-mpich=dir],[Specify MPICH install directory])], , with_mpich=none)
  if test "$with_mpich" != "none" ; then
    MPICH_DIR=$with_mpich;
    MPICH_CFLAGS=-I$MPICH_DIR/include
    MPICH_LDFLAGS=-L$MPICH_DIR/lib
    enable_mpi="mpich"
  fi

# OpenMPI // mpichが指定されておらず，かつ，ompが指定されていれば変数をセット
  if test "$with_mpich" = "none" ; then
    AC_ARG_WITH(ompi, [AC_HELP_STRING([--with-ompi=dir],[Specify OpenMPI install directory])], , with_ompi=none)

    if test "$with_ompi" != "none" ; then
      MPICH_DIR=$with_ompi;
      MPICH_CFLAGS=-I$MPICH_DIR/include
      MPICH_LDFLAGS=-L$MPICH_DIR/lib
      enable_mpi="ompi"
    fi
  fi

fi



#
# MPI Library for fortran
#

if test "$with_comp" = "FJ" ; then
  MPICH_LIBS=""
  case "$host" in
  sparc*)
    MPICH_FCLIBS="$MPICH_LIBS"" -lmpi_f77 -lmpi_f90"
    ;;
  esac
elif test "$enable_mpi" = "mpich" ; then
  MPICH_LIBS="-lmpich"
  MPICH_FCLIBS="$MPICH_LIBS"" -lfmpich"
  case "$TGT_DEVICE" in
  "Lion"|"Snow_Leopard"|"Leopard"|"Tiger_Intel"|"Tiger_PPC"|"Panther"|"Mountain_Lion")
#    MPICH_LIBS="$MPICH_LIBS"" -lpmpich"
    ;;
  esac
elif test "$enable_mpi" = "ompi" ; then
  MPICH_LIBS="-lmpi -lmpi_cxx"
  MPICH_FCLIBS="$MPICH_LIBS"" -lmpi_f77 -lmpi_f90"
  case "$TGT_DEVICE" in
  "Lion"|"Snow_Leopard"|"Leopard"|"Tiger_Intel"|"Tiger_PPC"|"Panther"|"Mountain_Lion")
#    MPICH_LIBS="$MPICH_LIBS"" -lmpi_cxx"
    ;;
  esac
fi


# warning. if a wrapper compiler is used, "enable_mpi" is empty, but parallel.
#
if test "$enable_mpi" = "none" ; then
  echo "  Warning: Neither MPICH nor OpenMPI library is not specified."
  echo "           If compilation error occurs about MPI, specify --with-mpich or --with-ompi option."
fi



#
# f90 example
#
AC_ARG_WITH(f90example, [AC_HELP_STRING([--with-f90example=(yes|no)],[make f90 example [yes]])], , with_f90example=yes)


#
# f90 compiler check
#
if test "_$FC" = "_"; then
  if test "_$F90" = "_"; then
    with_f90example="no"
  else
    FC=$F90
  fi
fi
if test "_$FCFLAGS" = "_"; then
  if test "_$F90FLAGS" != "_"; then
    FCFLAGS=$F90FLAGS
  fi
fi



#
# TextParser environment
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)
AC_SUBST(TP_LIBS)
if test "$with_parser" != "none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=`$TP_DIR/bin/tp-config --cflags`
  TP_LDFLAGS=`$TP_DIR/bin/tp-config --libs`
else
  echo "  Error: TextParser library is not specified."
  echo "         Specify --with-parser option."
  exit 1
fi



#
# PMlib environment
#
AC_SUBST(PM_DIR)
AC_SUBST(PM_CFLAGS)
AC_ARG_WITH(pm, [AC_HELP_STRING([--with-pm=dir],[Specify PMlib install directory])], , with_pm=none)
AC_SUBST(PM_LDFLAGS)

if test "$with_pm" != "none" ; then
  PM_DIR=$with_pm;
  PM_CFLAGS=`$PM_DIR/bin/pm-config --cflags`
  PM_LDFLAGS=`$PM_DIR/bin/pm-config --libs`
else
  echo "  Error: PM library is not specified."
  echo "         Specify --with-pm option."
  exit 1
fi




#
# CPM special flags
#
CPM_CFLAGS="-I$CPM_INST_DIR/include"
AC_SUBST(CPM_CFLAGS)

CPM_LDFLAGS="-L$CPM_INST_DIR/lib"
AC_SUBST(CPM_LDFLAGS)

CPM_LIBS="-lCPM"
AC_SUBST(CPM_LIBS)

CPM_CXX="$CXX"
AC_SUBST(CPM_CXX)

CPM_FC="$FC"
AC_SUBST(CPM_FC)

CPM_FC_LD="$FC"
AC_SUBST(CPM_FC_LD)

if test "$with_comp" = "FJ"; then
  CPM_FC_LD="$CXX"" --linkfortran"
elif test "$with_comp" = "INTEL" ; then
  CPM_LIBS="$CPM_LIBS"" -lstdc++ -lifport -lifcore"
fi
#CPM_LIBS="$CPM_LIBS"" -lstdc++"


#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc Examples/cxx"
if test "$with_f90example" = "yes" ; then
  MAKE_SUB_DIRS="$MAKE_SUB_DIRS"" Examples/f90"
fi


#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB


#
# Checks for header files.
#
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h \
                  string.h \
                  sys/time.h\
                  unistd.h \
                 ])


#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME



#
# Checks for library functions.
#
AC_CHECK_FUNCS([gethostname \
                gettimeofday \
                memset \
                strdup \
               ])


AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 doc/Makefile \
                 cpm-config \
                 include/cpm_Version.h \
                 Examples/cxx/Makefile \
                 Examples/f90/Makefile \
                ])

AC_OUTPUT

chmod +x ./cpm-config

