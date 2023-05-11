#!/bin/sh

#eval "$(conda shell.bash hook)"
#conda activate /eos/project/c/cms-ecal-calibration/ecal-conda/h4analysis
source /cvmfs/sft.cern.ch/lcg/releases/LCG_100/ROOT/v6.24.00/x86_64-centos7-gcc8-opt/ROOT-env.sh

export LD_LIBRARY_PATH=./lib:DynamicTTree/lib/:CfgManager/lib/:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=./interface:DynamicTTree/interface/:CfgManager/interface/:$ROOT_INCLUDE_PATH
