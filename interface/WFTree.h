#ifndef __WF_TREE__
#define __WF_TREE__

#include <string>
#include <vector>

#include "TTree.h"
#include "TString.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class WFTree
{
public: 
    //---ctors---
    WFTree(){};
    WFTree(int nSamples, uint64* idx, TTree* tree=NULL, string suffix="");
    //---dtor---
    ~WFTree(){};

    //---utils---
    void Init();
    void Fill();

    TTree* tree_; 
    string suffix_;

    uint64*       index;
    int           WF_samples;
    vector<int>   WF_ch; 
    vector<float> WF_time;
    vector<float> WF_val;
};

#endif
