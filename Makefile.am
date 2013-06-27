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

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS=cpm-config
EXTRA_DIST = cpm-config.in cpm-uname config_cpm.sh src/Makefile_spec LICENSE \
             doc/doxygen/Doxyfile doc/doxygen/makepdf.sh
dist_data_DATA = AUTHORS COPYING ChangeLog README LICENSE NEWS

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

