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

/**
   The TrackTree stores the information of all reconstructed tracks. It is closely
   related to the various PositionTree since the latters store the inputs of the 
   track fitting.
 */
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
    
    /**
       Number of tracks reconstructed in the event.
     */
    int                       n_tracks;
    /**
       X coordinates of each track.
     */
    std::vector<float>        X;
    /**
       Y coordinates of each track.
     */
    std::vector<float>        Y;    
    /**
       Number of hits of each track.
     */
    std::vector<int>          trackHits;
    /**
       Track fitting Chi2
     */
    std::vector<float>        trackChi2;
    /**
       Tracker layers that provided hits for this track.
     */
    std::vector<unsigned int> trackPattern;
    /**
       Minuit CovMatrixStatus.
     */
    std::vector<int>          fitStatus;
    /**
       Full TrackPar for each track.
     */
    std::vector<TrackPar>     fitResult;
};

#endif
