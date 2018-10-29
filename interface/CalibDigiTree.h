#ifndef __TRACK_TREE__
#define __TRACK_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;

class CalibDigiTree
{
public:


    //---ctors---
    CalibDigiTree(){};
    CalibDigiTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~CalibDigiTree(){};
  
    //---utils---
    void Init();
    void Clear() 
        {
	  n_calibWf=0;
	  wfChi2.clear();
        };

    void Fill() {tree_->Fill();};
  
    TTree*  tree_; 
    uint64* index;

    int                       n_calibWf;    
    std::vector<float>        wfChi2;
};

#endif
