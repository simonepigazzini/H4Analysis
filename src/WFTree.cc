#include "interface/WFTree.h"

WFTree::WFTree(int nSamples, uint64* idx, TTree* tree, string suffix)
{
    suffix_= suffix;
    tree_ = tree ? tree : new TTree();

    index=idx;
    WF_samples = nSamples;
}

void WFTree::Init()
{
    //---set total number of WF samples 
    WF_ch.reserve(WF_samples);
    WF_time.reserve(WF_samples);
    WF_val.reserve(WF_samples);
    //---global branches
    tree_->Branch("index", index,"index/l");
    tree_->Branch(("WF_ch"+suffix_).c_str(), &WF_ch);
    tree_->Branch(("WF_time"+suffix_).c_str(), &WF_time);
    tree_->Branch(("WF_val"+suffix_).c_str(), &WF_val);
}

void WFTree::Fill()
{
    //---fill tree
    tree_->Fill();

    //---automatically reset after filling
    Reset();
}

void WFTree::Reset()
{
    //---clear containers but preserve size
    WF_ch.clear();
    WF_time.clear();
    WF_val.clear();
    WF_ch.reserve(WF_samples);
    WF_time.reserve(WF_samples);
    WF_val.reserve(WF_samples);
}        
