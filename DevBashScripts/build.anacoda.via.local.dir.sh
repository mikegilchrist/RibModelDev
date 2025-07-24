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

if [[ "$#" -eq 0 ]]; then
    echo "sourceDir unspecified, using default value: \"./RibModelFramework\"";
    sourceDir="RibModelFramework";
else
    echo "Using sourceDir: $1"
    sourceDir="$1"
fi

rVersion=$(Rscript -e 'cat(paste0(R.version$major, ".", sub("\\.[0-9]+", "", R.version$minor)))')
echo "R Version: $rVersion"
branch=$(git -C "$sourceDir" branch --show-current)
echo "Git Branch: $branch"
#installDir="~/R/lib/$rVersion-dev/"
installDir=~/R/lib/"$branch"
# ensure directory exists
mkdir -p "$installDir"
echo "Install Dir: $installDir"




if R CMD build --no-build-vignettes "$sourceDir"; then
    echo "Build succeeded. Installing package";
    MAKE="make -j$(($(nproc)-1))"; #use 1 less than number of cores on machine
    if R CMD INSTALL -l "$installDir" $(ls -t AnaCoDa_*.tar.gz| head -1) ; then
	echo "Install succeeded."
	## test code
	#cd "$sourceDir/tests";
	## keep this one long command (note the double quotes
	if Rscript -e ".libPaths(\"$installDir\"); 
	   .libPaths();
           setwd(\"$sourceDir/tests/\");
	   source(\"./testthat.R\");
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
