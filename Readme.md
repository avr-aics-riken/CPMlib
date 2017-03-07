# CPMlib - Computational space Partitioning Management library

* Copyright (c) 2012-2014 Institute of Industrial Science (IIS), The University of Tokyo. All rights reserved.
* Copyright (c) 2014-2016 Advanced Institute for Computational Science (AICS), RIKEN. All rights reserved.
* Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University. All rights reserved.

## OUTLINE

CPMlib is a C++ class library to efficiently assist the development of unsteady physical simulators. This library provides functions, such as the management of partitioned subdomains by domain decomposition method and wrapper functions for neighbor and global communications.

## SOFTWARE REQUIREMENT
- Cmake
- MPI library
- TextParser
- PMlib(option for a code array example)


## INGREDIENTS
~~~
ChangeLog.md      History of development
License.txt       License to apply
Readme.md         This document, including the description of build
cmake/            Modules of cmake
doc/              Documents
Examples/         Example source codes
include/          Header files
src/              Source codes
~~~

## HOW TO BUILD

### Build

~~~
$ export PM_HOME=/hogehoge
$ mkdir BUILD
$ cd BUILD
$ cmake [options] ..
$ make
$ sudo make install
~~~


### Options

`-D INSTALL_DIR=install_directory`

>  Specify the directory that this library will be installed. Built library is installed at `install_directory/lib` and the header files are placed at `install_directory/include`. The default install directory is `/usr/local/PMlib`.

`-D real_type=` {float | double}

>  Specify the type of floating point. If this option is omitted, the default is float.

`-D enable_LMR=` {no | yes}

>  This option allows you to use LMR framework. The default is no.

`-D with_MPI=` {yes | no}

>  If you use an MPI library, specify `with_MPI=yes`, the default is yes.

`-D with_example=` {no | yes}

>  This option turns on compiling sample codes. The default is no.

`-D with_TP=` *Installed_directory*

> Specify the directory path that TextParser is installed.




The default compiler options are described in `cmake/CompilerOptionSelector.cmake` file. See BUILD OPTION section in CMakeLists.txt in detail.


## Configure Examples

`$ export CMP_HOME=hogehoge`

In following exsmples, assuming that TextParser and PMlib are installed under the CPM_HOME directory. If not, please specify applicable directory paths.

### INTEL/GNU compiler

~~~
$ cmake -DINSTALL_DIR=${CPM_HOME}/CPMlib -Dwith_MPI=yes -Dreal_type=float -Denable_LMR=no -Dwith_example=no -Dwith_TP=${CPM_HOME}/TextParser ..
~~~


### FUJITSU compiler / FX10 on login nodes (Cross compilation)

~~~
$ cmake -DINSTALL_DIR=${CPM_HOME}/CPMlib -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx10.cmake -Dwith_MPI=yes -Dreal_type=float -Denable_LMR=no -Dwith_example=no -Dwith_TP=${CPM_HOME}/TextParser ..
~~~


### FUJITSU compiler / K computer on login nodes (Cross compilation)

~~~
$ cmake -DINSTALL_DIR=${CPM_HOME}/CPMlib -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake -Dwith_MPI=yes -Dreal_type=float -Denable_LMR=no -Dwith_example=no -Dwith_TP=${CPM_HOME}/TextParser ..
~~~


##### Note
- Before building, execute following command for clean. `$ make distclean`


### Other platforms

Refer to `Makefile_hand` hand for manually prepare the Makefiles for
a system not covered in the previous sections.
Edit environmental variables in a `Makefile_hand` file for a target machine.


## EXAMPLES

* If you specify the test option by `-Denable_example=yes`, you can
execute the intrinsic tests by;

	`$ make test` or `$ ctest`

* The detailed results are written in `build/Testing/Temporary/LastTest.log` file.
Meanwhile, the summary is displayed for stdout.




## CONTRIBUTORS

* Kenji     Ono      *keno@{cc.kyushu-u.ac, riken}.jp*
* Yasuhiro  Kawashima
