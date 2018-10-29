#include "interface/CalibDigiTree.h"

CalibDigiTree::CalibDigiTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();
    n_calibWf = 0;
    index=idx;
}

void CalibDigiTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");
    tree_->Branch("n_calibWf", &n_calibWf, "n_calibWf/I");
    tree_->Branch("wfChi2", &wfChi2);
}
