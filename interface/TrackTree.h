#ifndef __TRACK_TREE__
#define __TRACK_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;

//---Version of track class for persistent storage in the output file
class TrackPar
{
public:
    TrackPar() {}
    ~TrackPar() {}

    double x() { return  value[0]; }
    double y() { return  value[1]; }
    double alpha() { return  value[2]; }
    double beta() { return  value[3]; }

    double err_x() { return  sqrt(covariance[0]); }
    double err_y() { return  sqrt(covariance[2]); }
    double err_alpha() { return  sqrt(covariance[5]); }
    double err_beta() { return  sqrt(covariance[9]); }

    double corr_x_alpha()  { return covariance[3]/(err_x() * err_alpha()); }
    double corr_y_beta()  { return covariance[7]/(err_y() * err_beta()); }
    double corr_x_y() { return covariance[1]/(err_x() * err_y()); }
    double corr_alpha_beta() { return covariance[8]/(err_alpha() * err_beta()); }

    std::vector<double> value;
    std::vector<double> covariance;

//    ClassDef(TrackPar, 1)
};

class TrackTree
{
public:


    //---ctors---
    TrackTree(){};
    TrackTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~TrackTree(){};
  
    //---utils---
    void Init();
    void Clear() 
        {
            X.clear();
            Y.clear();            
            trackHits.clear();
            trackPattern.clear();
            trackChi2.clear();
            fitStatus.clear();
            fitResult.clear();
        };

    void Fill() {tree_->Fill();};
  
    TTree*  tree_; 
    uint64* index;
    
    int                       n_tracks;
    std::vector<float>        X;
    std::vector<float>        Y;    
    std::vector<int>          trackHits;
    std::vector<float>        trackChi2;
    std::vector<unsigned int> trackPattern;
    std::vector<int>          fitStatus;
    std::vector<TrackPar>     fitResult;
};

#endif
