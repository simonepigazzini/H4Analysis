#include "interface/PositionTree.h"

PositionTree::PositionTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();
    
    index=idx;
}

void PositionTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");

    //---position branches
    n_clusters = 0;
    clusters = vector<PositionMeasurement>();
    tree_->Branch("n_clusters", &n_clusters, "n_clusters/I");
    tree_->Branch("clusters", &clusters);
}

void PositionTree::Clear()
{
    n_clusters = 0;
    clusters.clear();
}
