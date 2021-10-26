#ifndef __POSITION_TREE__
#define __POSITION_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;


//---General position measurement used for permanent storage in output file.  
class PositionMeasurement
{
public:
    //---ctors
    PositionMeasurement():
        nHits_(0),
        X_(0),
        Y_(0),
        magnitude_(-1)
        {};
    PositionMeasurement(int nhits, float x, float y, float mag=0):
        nHits_(nhits),
        X_(x),
        Y_(y),
        magnitude_(mag)
        {};
    
    /**
       Number of hits in the cluster.
     */
    int   nHits_;
    /**
       X coordinate of the cluster.
     */
    float X_;
    /**
       Y coordinate of the cluster.
     */
    float Y_;
    /**
       Cluster total amplitude.
     */
    float magnitude_;
};
    
/**
   The position tree contains all measurements from a specific tracking plane (in H4 one hodoscope
   plane). Each plane can have have multiple hits groupped in clusters. See PositionMeasurement
   for a description of the measurement.
 */
class PositionTree
{
public: 
    //---ctors---
    PositionTree(){};
    PositionTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~PositionTree(){};
    
    //---utils---
    void Init();
    void Clear();
    void Fill() {tree_->Fill();};    
    
    TTree*  tree_; 
    uint64* index;

    /**
       Number of clusters in a plane
     */
    int                         n_clusters;
    /**
       Clusters from this plane (see PositionMeasurement).
     */
    vector<PositionMeasurement> clusters;
};

#endif
