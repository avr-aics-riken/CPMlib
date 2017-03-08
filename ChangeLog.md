# CPMlib - Computaional space Partitioning Management library

## REVISION HISTORY

---
- 2017-03-07 Version 2.3.3
  - modify implementation of 'real_type'


---
- 2017-03-07 Version 2.3.2
  - add fx100 entry
  - In case of Intel compiler, change linker  from CXX to Fortran for f90 example
  - Tested.

  |Compiler|Serial|Ex.|MPI |Ex.|
  |:--|:--:|:--:|:--:|:--:|
  |Intel 17.0.1 |||ok|100%|
  |GNU 6.2.0    |||ok|66%|
  |fx10         |||ok|–|
  |fx100        |||ok|–|
  |K            |||ok|–|


---
- 2017-2-23 Version 2.3.1
  - remove direction that relates PMlib for a CPMlib.


---
- 2017-2-23 Version 2.3.0
  - modify Toolchain_K.cmake so that user can build with cmake version 2.6 on K.
  - suppress OpenMP option building for a CPMlib itself, but still enable for examples.


---
- 2017-2-16 Version 2.2.4
  - enable CPM_VERSION


---
- 2017-2-15 Version 2.2.3
  - bug fix : error occurs when GNU compiler is used.
  	 - `src/cpm_ParaManager_frtIF.cpp`, remove arg `false` at lines 2066, 3138, 3277
  - If users employ openmpi-2.x.x, configure option `--enable-mpi-cxx` must be given to enable C++ binding.


---
- 2017-2-12 Version 2.2.2
  - modify Readme.md
  - correct fx10 compiler options


---
- 2017-2-11 Version 2.2.1
  - modify Readme.md


---
- 2017-2-7 Version 2.2.0
  - cmake branch in github
  - change acronym
  - LMR option
  - introduce getarg() for mconvp_CPM.f90


---
- 2016-10-31 Version 2.1.5
  - add LMR function and padding


---
- 2015-11-28 Version 2.1.2
	- r38 modify LDFLAGS option of GNU compiler in INSTALL


---
- 2015-11-28 Version 2.1.1
	- r37 --with-comp becomes essential
	- change configure.ac, INSTALL, and NEWS


---
- 2015-07-25 Version 2.1.0
	- r36 Intel MPI option
	- add --with-impi option


---
- 2015-07-08 Version 2.0.5
	- r35 bug fix when building by gnu compiler


---
- 2015-06-22 Version 2.0.4
	- r34 add Global2LocalIndex()

---
- 2015-06-15 Version 2.0.3
	- r33 modify examle options in configure


---
- 2015-06-10 Version 2.0.2
	- r32 clean package


---
- 2015-06-09 Version 2.0.1
	- r31 introduce BUILD_DIR to keep source directory clean
	- Change to run configure
	- Change configure.ac


---
- 2015-06-05 Version 2.0.0
	- r30 LMR branch
	- version 2.0.0


---
- 2015-03-14 Version 1.2.4
	- r29 update year
	- add mpiicc, mpiicpc
	- update description of INSTALL, add `${FFV_HOME}`


---
- 2014-12-25 Version 1.2.3   


---
- 2014-03-19
	- r28 update user guide


---
- 2014-03-19 Version 1.2.2
	- r27 suppress f90 example [ver 1.2.2]
	- add Examples in `EXTRA_DIST`


---
- 2014-03-19 Version 1.2.1
	- r26 enhance algorithm for domain division [ver 1.2.1]
	- add enum cpm_DivPolicy
	- divide `DecideDivPattern()` to `DecideDivPattern_CommSize()` and `DecideDivPattern_Cube()`
	- change `cpm_ParaManager::VoxelInit() (add argument [divPolicy])`
  	- add AICS copyright


---
- 2014-03-13 Version 1.2.0
	- r25 bug fix for ver 1.1.9 [ver 1.2.0]


---
- 2014-03-12 Version 1.1.9
	- r24 enhance algorithm for domain division [ver 1.1.9]


---
- 2014-03-12 Version 1.1.8
	- r23 update domain division algorithm [ver 1.1.8]
	- add `CPM_ERROR_DECIDE_DIV_PATTERN`
	- change `CalcCommSize() iDiv >> iDiv-1`
	- change `DecideDivPattern()`


---
- 2014-03-04 Version 1.1.7
	- r22 version format [ver 1.1.7]
	- change output format for cpm-config --version


---
- 2014-02-11 Version 1.1.6
	- r21 update configure.ac [ver 1.1.6]
	- remove chk-uname, config-cpm.sh
	- add GNU compiler
	- add compiler environment CC for NERSC Hopper


---
- 2013-10-30 Version 1.1.5
	- r20 update chk-uname [ver 1.1.5]


---
- 2013-10-02 Version 1.1.4
	- r19 modify configure.ac for CPM special flag [ver 1.1.4]
	- `CPM_LIBS="$CPM_LIBS"" -lstdc++"`


---
- 2013-10-02 Version 1.1.3
	- r18 modify for intel mpi [ver 1.1.3]
	- include mpi.h before stdio.h to suppress error message `#error "SEEK_SET is #defined but must not be for the C++ binding of MPI"`


---
- 2013-07-20 Version 1.1.2
	- r17 void VersionInfo() >> std::string getVersionInfo() [ver 1.1.2]
	- Also change Examples/cxx/main()


---
- 2013-07-13
	- r16 update `include/cpm_Version.h.in` and cpm_Base.h
	- cpm_Version.h.in to generate `CPM_REVISION` automatically
	- change Versioninfo()


---
- 2013-07-13 Version 1.1.1
	- r15 update configure.ac [ver 1.1.1]
	- modify link for compiling with OpenMPI 1.6 or later, which requires -lmpi_cxx
	- modify with-comp option if FJ else >> if FJ elif INTEL


---
- 2013-07-12 Version 1.1.0
	- r14 update INSTALL
	- modify example of configure for K


---
- 2013-07-06
	- r13 update configure [ver 1.1.0]
	- cpm-uname -> chk-uname
	- update ChangeLog, NEWS
	- change behavior when a wrapper compiler is specified
	- remove extra bused setting in configure.ac


---
- 2013-06-29 Version 1.0.12
	- r12 exp. branch, introduce usr_TPdomain()
	- Text implementation of pure virtual method for division policy (experimental branch）
	- r11 update to 1.0.12


---
- 2013-06-27 Version 1.0.11
  - r10 ver. 1.0.11


---
- 2013-06-27 Version 1.0.10
  - r9 v1.0.10


------------------
- 2013-06-25 Version
  - r8 update README
  - r7 update AUTHORS


------------------
- 2013-06-22
   - r6 improve 1.0.9


------------------
- 2013-06-19
  - r5 CPMlib-1.0.9


---
- 2013-05-02 Version 1.0.9


---
- 2013-04-03 Version 1.0.8   


---
- 2012-08-25 Version 1.0.6


---
- 2012-08-04 Version 1.0.5


---
- 2012-07-09 Version 1.0.4   
  - r4 CPMlib-1.0.8


---
- 2012-07-03 Version 1.0.3
  - r3 CPMlib-1.0.6


---
- 2012-06-27 Version 1.0.2
  - r2 CPMlib-1.0.2


---
- 2012-06-25 Version 1.0.1   


---
- 2012-06-18 Version 1.0.0
