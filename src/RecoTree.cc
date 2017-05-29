#include "interface/RecoTree.h"

RecoTree::RecoTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();

    index=idx;
    start_time=0;
    time_stamp=0;
    evt_flag=0;
    run=0;
    spill=0;
    event=0;   
    //---global branches
    tree_->Branch("index", index, "index/l");
    tree_->Branch("start_time", &start_time, "start_time/l");
    tree_->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree_->Branch("evt_flag", &evt_flag, "evt_flag/i");    
    tree_->Branch("run", &run, "run/i");
    tree_->Branch("spill", &spill, "spill/i");
    tree_->Branch("event", &event, "event/i");
}

RecoTree::~RecoTree()
{
    delete tree_;
}
