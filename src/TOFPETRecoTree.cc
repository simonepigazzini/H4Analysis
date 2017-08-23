#include "interface/TOFPETRecoTree.h"

TOFPETRecoTree::TOFPETRecoTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();
    
    index=idx;
}

void TOFPETRecoTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");
    //---tofpet branches
    tree_->Branch("t_sipm", &t_sipm, "t_sipm/D");
    tree_->Branch("tot", &tot, "tot/D");
    tree_->Branch("energy", &energy, "energy/D");        
    tree_->Branch("t_h4daq", &t_h4daq, "t_h4daq/l");
    tree_->Branch("t_tofpet", &t_tofpet, "t_tofpet/l");
}
