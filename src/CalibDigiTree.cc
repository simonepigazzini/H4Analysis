#include "interface/CalibDigiTree.h"

CalibDigiTree::CalibDigiTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();
    n_calibWf = 0;
    index=idx;
}

void CalibDigiTree::Init()
{
    board = -1;
    group = -1;
    channel = -1;

    calib = new DigiChannelCalibration();
    
    //---global branches
    tree_->Branch("index", index, "index/l");
    tree_->Branch("n_calibWf", &n_calibWf, "n_calibWf/I");
    tree_->Branch("wfChi2", &wfChi2);
    tree_->Branch("board", &board, "group/I");
    tree_->Branch("group", &group, "board/I");
    tree_->Branch("channel", &channel, "channel/I");
    tree_->Branch("calib", calib);

    return;
}
