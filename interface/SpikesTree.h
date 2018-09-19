#ifndef __SPIKES_TREE__
#define __SPIKES_TREE__

#include <memory>
#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class SpikesTree
{
public: 
    //---ctors---
    SpikesTree(){};
    SpikesTree(uint64* idx, TTree* tree=NULL, string prefix="");
    //---dtor---
    ~SpikesTree();

    //---utils---
    void Init(vector<string>& names);
    void Fill() {tree_->Fill();};
    
    TTree*        tree_; 
    string        prefix_;

    uint64*       index;
    unsigned int  n_channels;
    int           max_hit;
    int*          channels;
    float*        undershoot;
    float*        amp_sum_matrix;    
};

#endif
