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

rVersion="3.6"
installDir="~/R/x86_64-pc-linux-gnu-library/$rVersion"
if R CMD build RibModelFramework ; then
    echo "Build succeeded. Installing package";
    MAKE="make -j$(($(nproc)-1))"; #use 1 less than number of cores on machine
    if R CMD INSTALL -l $installDir $(ls -t AnaCoDa_*.tar.gz| head -1) ; then
	echo "Install succeeded."
	## test code
	cd RibModelFramework/tests;
	## keep this one long command (note the double quotes
	if Rscript -e ".libPaths(\"$installDir\"); 
	   .libPaths(); \\
	   source(\"testthat.R\"); \\
	   packageVersion(\"AnaCoDa\")"; 
	then
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
