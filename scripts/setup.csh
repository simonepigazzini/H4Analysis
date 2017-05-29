#!/bin/csh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.csh
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/6.05.02/x86_64-slc6-gcc49-opt/root/
source ./bin/thisroot.csh
cd -

setenv LD_LIBRARY_PATH DynamicTTree/lib/:$LD_LIBRARY_PATH
