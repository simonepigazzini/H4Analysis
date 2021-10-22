#ifndef __TRG_TREE__
#define __TRG_TREE__

#include <string>
#include <vector>

#include "TTree.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class TrgTree
{
public: 
    //---ctors---
    TrgTree(){};
    TrgTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~TrgTree(){};

    //---utils---
    void Init(map<int, string>& triggers);
    void Fill() {tree_->Fill();};
    void Reset();
    
    TTree* tree_; 

    uint64*       index;
    vector<int>   trgs_;
    unsigned int  trg_type_;
};

#endif
