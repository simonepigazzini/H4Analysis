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
    n_clusters_X = 0;
    n_clusters_Y = 0;
    cluster_X_size = vector<int>(); 
    cluster_Y_size = vector<int>(); 
    X = vector<float>(); 
    Y = vector<float>(); 
    tree_->Branch("n_clusters_X", &n_clusters_X, "n_clusters_X/I");
    tree_->Branch("n_clusters_Y", &n_clusters_Y, "n_clusters_Y/I");    
    tree_->Branch("cluster_X_size", &cluster_X_size);
    tree_->Branch("cluster_Y_size", &cluster_Y_size);
    tree_->Branch("X", &X);
    tree_->Branch("Y", &Y);
}
