H4Analysis
==========
![docs](https://gitlab.cern.ch/spigazzi/H4Analysis/badges/master/pipeline.svg)

# [Documentation](https://h4analysis.web.cern.ch/)
  - This repository aims to provide a fast reconstruction of data
    acquired with the H4DAQ.
  - The main executable is =bin/H4Reco=. The data processing is 
    steered through the input cfg file (examples in cfg/).

# Requirements
  - c++17
  - ROOT

# Install and run
   git clone --recursive https://github.com/simonepigazzini/H4Analysis.git

   cd H4Analysis

   source script/setup.sh(csh) (only on lxplus, or any machine with afs and EOS access).

   make

   bin/H4Reco cfg/Oct2015_timing.cfg 4443
     
