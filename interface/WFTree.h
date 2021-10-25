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

/**
   Pulse shapes (waveforms) tree:
   - This tree holds a copy of the pulse shapes stored in the digi branches of the raw 
     data tree.
   - Usually the entries of this tree are prescaled (see WFAnalyzer).
   - each sample has a time (#WF_time) and a amplitude (#WF_val).
   - Samples from all channels are stored in the same branches, use #WF_ch to filter out
     samples from a specific channel.
 */
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
    void Reset();
    
    TTree* tree_; 
    string suffix_;

    uint64*       index;
    int           WF_samples;    
    /**
       Sample channel (channel name<->index map is stored in the DigiTree)
     */
    vector<int>   WF_ch; 
    /**
       Sample time (channel name<->index map is stored in the DigiTree)
     */
    vector<float> WF_time;
    /**
       Sample amplitude (channel name<->index map is stored in the DigiTree)
     */
    vector<float> WF_val;
};

#endif
