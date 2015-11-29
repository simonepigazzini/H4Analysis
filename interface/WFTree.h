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
    WFTree(int nCh, int nSamples, int* idx, TString suffix="", Long64_t autosave=10);
    //---dtor---
    ~WFTree();

    //---utils---
    void Fill() {tree_->Fill();};
    void Write(const char* name="wf_tree", const char* title="wf_tree")
        {tree_->BuildIndex("index"); tree_->SetTitle(title); tree_->Write(name);};

    TTree* tree_; 
    TString suffix_;

    int*   index;
    uint64 time_stamp;
    int    WF_samples;
    int*   WF_ch; 
    float* WF_time;
    float* WF_val;
};

WFTree::WFTree(int nCh, int nSamples, int* idx, TString suffix, Long64_t autosave)
{
    suffix_= suffix;
    tree_ = new TTree("wf", "wf_tree");
    tree_->SetMaxVirtualSize(10000000);

    index=idx;
    time_stamp=0;
    //---set total number of WF samples 
    WF_samples = nSamples*nCh;
    WF_ch = new int[WF_samples];
    WF_time = new float[WF_samples];
    WF_val = new float[WF_samples];
    //---global branches
    tree_->Branch("index", index,"index/I");
    tree_->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree_->Branch(Form("WF_samples%s",suffix_.Data()), &WF_samples, Form("WF_samples%s/I",suffix_.Data()));
    tree_->Branch(Form("WF_ch%s",suffix_.Data()), WF_ch,            Form("WF_ch%s[WF_samples%s]/I",suffix_.Data(),suffix_.Data()));
    tree_->Branch(Form("WF_time%s",suffix_.Data()), WF_time,        Form("WF_time%s[WF_samples%s]/F",suffix_.Data(),suffix_.Data()));
    tree_->Branch(Form("WF_val%s",suffix_.Data()), WF_val,          Form("WF_val%s[WF_samples%s]/F",suffix_.Data(),suffix_.Data()));
}

WFTree::~WFTree()
{
    delete[] WF_ch; 
    delete[] WF_time;
    delete[] WF_val;

    delete tree_;
}

#endif
