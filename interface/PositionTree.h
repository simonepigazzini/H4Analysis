#ifndef __POSITION_TREE__
#define __POSITION_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;

class PositionTree
{
public: 
    //---ctors---
    PositionTree(){};
    PositionTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~PositionTree(){};
    
    //---utils---
    void Init();
    void Fill() {tree_->Fill();};

    TTree*  tree_; 

    uint64*        index;
    int            n_clusters_X;
    int            n_clusters_Y;    
    vector<int>    cluster_X_size;
    vector<int>    cluster_Y_size;
    vector<float>  X;
    vector<float>  Y;
};

#endif
