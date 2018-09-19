#ifndef __FITPIX_TREE__
#define __FITPIX_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;

/* #define FITPIX_MAX_HITS 5000 */
/* #define FITPIX_MAX_CLUSTERS 5000 */

class FitpixTree
{
public: 
    //---ctors---
    FitpixTree(){};
    FitpixTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~FitpixTree(){};
    
    //---utils---
    void Init();
    void Clear() 
    {
      hitX.clear();
      hitY.clear();
      hitCharge.clear();

      clusterX.clear();
      clusterY.clear();
      clusterCharge.clear();
      clusterSize.clear();
    };

    void Fill() {tree_->Fill();};

    TTree*  tree_; 

    uint64* index;
    int n_hits;
    /* float* hitX; */
    /* float* hitY; */
    /* float* hitCharge; */
    std::vector<float> hitX;
    std::vector<float> hitY;
    std::vector<float> hitCharge;
    
    int n_clusters;
    /* float* clusterX; */
    /* float* clusterY; */
    /* float* clusterCharge; */
    /* int* clusterSize; */
    std::vector<float> clusterX;
    std::vector<float> clusterY;
    std::vector<float> clusterCharge;
    std::vector<int> clusterSize;    
};

#endif
