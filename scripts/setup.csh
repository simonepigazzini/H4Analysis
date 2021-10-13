#!/bin/csh

conda activate /eos/project/c/cms-ecal-calibration/ecal-conda/h4analysis
    
setenv LD_LIBRARY_PATH DynamicTTree/lib/:$LD_LIBRARY_PATH
