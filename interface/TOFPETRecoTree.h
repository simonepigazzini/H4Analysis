#ifndef __TOFPET_RECO_TREE__
#define __TOFPET_RECO_TREE__

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;

class TOFPETRecoTree
{
public: 
    //---ctors---
    TOFPETRecoTree(){};
    TOFPETRecoTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~TOFPETRecoTree(){};
    
    //---utils---
    void Init();
    void Fill() {tree_->Fill();};

    TTree*  tree_; 

    uint64* index;
    double t_sipm;
    double tot;
    double energy;
    uint64 t_h4daq;
    uint64 t_tofpet;
};

#endif

    
