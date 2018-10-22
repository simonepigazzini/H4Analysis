#ifndef __POSITION_TREE__
#define __POSITION_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;


//---General position measurement used for permanent storage in output file  
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
    
    inline int   nHits() {return nHits_;};
    inline float X() {return X_;};
    inline float Y() {return Y_;};
    inline float magnidute() {return magnitude_;};

private:
    int   nHits_;
    float X_;
    float Y_;
    float magnitude_;
};
    
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
    
    int                         n_clusters;
    vector<PositionMeasurement> clusters;
};

#endif
