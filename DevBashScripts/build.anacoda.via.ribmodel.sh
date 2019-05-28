#!/bin/bash

# Script to build and install AnaCoDa package from a git repository
# Should be run above git repo (I use symbolic links in my repo directory)

#You may need to reinstall Rcpp
#If you want a subwindow in figures you need to install packages
##  survival
##  Matrix
##  Hmisc

##build the package outside the repository RibModelFramework
## The tar.gz name for package is automatically set
#
# Global install
#R CMD build RibModelFramework; sudo MAKE="make -j4" R CMD INSTALL $(ls -t ribModel_*.tar.gz| head -1)
#
#  Local install

if R CMD build RibModelFramework ; then
    echo "Build succeeded. Installing package";
    MAKE="make -j8";
    if R CMD INSTALL -l ~/R/lib-dev $(ls -t AnaCoDa_*.tar.gz| head -1) ; then
	echo "Install succeeded."
	## test code
	cd RibModelFramework/tests;
	if Rscript -e ".libPaths(\"~/R/lib-dev\"); .libPaths(); source(\"testthat.R\"); packageVersion(\"AnaCoDa\")"; then
	    echo "TestThat succeeded"
	else
	    echo "TestThat failed"
	fi
    else
	echo "Install failed."
    fi
    
else
    echo "Build failed.  Try again"
fi
