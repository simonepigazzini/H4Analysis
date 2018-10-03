#include "interface/FitpixTree.h"

FitpixTree::FitpixTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();

    n_hits = 0;
    n_clusters = 0;
    index=idx;
}

void FitpixTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");

    // hitX = new float[FITPIX_MAX_HITS];
    // hitY = new float[FITPIX_MAX_HITS];
    // hitCharge = new float[FITPIX_MAX_HITS];

    // clusterX = new float[FITPIX_MAX_CLUSTERS];
    // clusterY = new float[FITPIX_MAX_CLUSTERS];
    // clusterCharge = new float[FITPIX_MAX_CLUSTERS];
    // clusterSize = new int[FITPIX_MAX_CLUSTERS];

    tree_->Branch("n_hits", &n_hits, "n_hits/I");
    tree_->Branch("n_clusters", &n_clusters, "n_clusters/I");
    // tree_->Branch("hitX", hitX ,"hitX[n_hits]/F");
    // tree_->Branch("hitY", hitY ,"hitY[n_hits]/F");
    // tree_->Branch("hitCharge", hitCharge ,"hitCharge[n_hits]/F");
    tree_->Branch("hitX", &hitX );
    tree_->Branch("hitY", &hitY );
    tree_->Branch("hitCharge", &hitCharge );

    // tree_->Branch("clusterX", clusterX, "clusterX[n_clusters]/F");
    // tree_->Branch("clusterY", clusterY, "clusterY[n_clusters]/F");
    // tree_->Branch("clusterCharge", clusterCharge, "clusterCharge[n_clusters]/F");
    // tree_->Branch("clusterSize", clusterSize, "clusterSize[n_clusters]/I");
    tree_->Branch("clusterX", &clusterX);
    tree_->Branch("clusterY", &clusterY);
    tree_->Branch("clusterCharge", &clusterCharge);
    tree_->Branch("clusterSize", &clusterSize);
}
