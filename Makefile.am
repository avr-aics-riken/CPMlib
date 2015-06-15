###
###################################################################################
#
# CPMlib - Cartesian Partition Manager Library
# 
# Copyright (C) 2012-2014 Institute of Industrial Science, The University of Tokyo.
# All right reserved.
#
# Copyright (c) 2014-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
###

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS=cpm-config
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = cpm-config.in src/Makefile_spec LICENSE BUILD_DIR/.gitkeep \
             doc/doxygen/Doxyfile doc/doxygen/makepdf.sh
dist_data_DATA = AUTHORS COPYING ChangeLog README LICENSE NEWS

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*
