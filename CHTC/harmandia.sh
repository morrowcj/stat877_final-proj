#!/bin/bash

# untar your R installation
tar -xzf R_ML2.tar.gz

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export RHOME=$(pwd)/R

# run R, with the name of your  R script
Rscript harmandia.r $1 

# pass seed 