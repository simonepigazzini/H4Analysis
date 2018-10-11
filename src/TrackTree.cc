#include "interface/TrackTree.h"

TrackTree::TrackTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();
    n_tracks = 0;
    index=idx;
}

void TrackTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");
    tree_->Branch("n_tracks", &n_tracks, "n_tracks/I");
    tree_->Branch("fitResult", &fitResult);
    tree_->Branch("fitStatus", &fitStatus);
    tree_->Branch("trackHits", &trackHits);
    tree_->Branch("trackChi2", &trackChi2);
}
