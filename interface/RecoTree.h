#ifndef __RECO_TREE__
#define __RECO_TREE__

#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

class RecoTree
{
public: 
    //---ctors---
    RecoTree(){};
    RecoTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~RecoTree();

    //---utils---
    void Fill() {tree_->Fill();};
    void Write(const char* name="h4", const char* title="main_h4reco_tree")
        {tree_->BuildIndex("index"); tree_->SetTitle(title); tree_->Write(name);};
    void AddFriend(const char* friend_name="wf_tree = wf") {tree_->AddFriend(friend_name, "");};

    TTree* tree_; 

    uint64*      index;
    uint64       start_time;
    uint64       time_stamp;
    unsigned int evt_flag;
    unsigned int run;
    unsigned int spill;
    unsigned int event;       
};

#endif
