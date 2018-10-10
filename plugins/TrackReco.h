#ifndef __TRACK_RECO__
#define __TRACK_RECO__

#include "interface/PluginBase.h"
//#include "interface/TrackTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "Math/SMatrix.h"
#include "Math/GenVector/Rotation3D.h"
using namespace std;

using MeasurementErrorMatrix_t = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >;
using RotationMatrix_t = ROOT::Math::Rotation3D;
using Measurement_t = TVector2;
using GlobalMeasurement_t = TVector3;
using GlobalCoord_t = TVector3;

class TrackReco: public PluginBase
{
public:

    //---ctors---
    TrackReco() {};
  
    //---dtor---
    ~TrackReco() {};

    class TrackLayer 
    {
    public:
    TrackLayer(const GlobalCoord_t& pos):
      position_(pos),rotation_()
      {
      };

    TrackLayer(const GlobalCoord_t& pos, const RotationMatrix_t& rot):
      position_(pos),rotation_(rot)
      {
      };
      
      inline void globalCoordinates( const Measurement_t& loc, GlobalCoord_t& pos ) const
      {
      };
      
      /*
	inline void localCoordinates( GlobalCoord_t& loc, GlobalCoord_t& pos ) const;  
	 { 
	 } 
      */
      
      GlobalCoord_t position_; //axis origin translation from localFrame to global Frame 
      RotationMatrix_t rotation_; //axis rotation matrix from localFrame to global Frame 
    };

    class TelescopeLayout
    {
    public:
      TelescopeLayout()
      {
	layers_.clear();
      };
      
      void addLayer(const TrackLayer& layer)
      {
	layers_.push_back(layer);
      }

      inline void globalCoordinates( const int& k, const Measurement_t& loc, GlobalCoord_t& pos ) const
      {
	layers_[k].globalCoordinates(loc, pos);
      }
            
      std::vector<TrackLayer> layers_;
    };
    
    class TrackMeasurement {
    public:
    TrackMeasurement( double x, double y, const TelescopeLayout& hodo, int layer ) :
      localPosition_(x,y), localPositionError_(), hodo_(hodo), layer_(layer)
      {
	
      }
      
      ~TrackMeasurement() {};
      
      inline void setVarianceX(double& x) {};
      inline void setVarianceY(double& y) {};

      inline void toGlobal( GlobalCoord_t& pos ) const
      {
	hodo_.globalCoordinates(layer_, localPosition_, pos);
      }
      
      
      Measurement_t localPosition_;
      MeasurementErrorMatrix_t localPositionError_;
      const TelescopeLayout& hodo_;
      int layer_;
    };
  
    class Track {
    public:
    Track(const TelescopeLayout& hodo) : hodo_(hodo)
      {
	hits_.clear();
      };
      
      ~Track() {};
      
      inline Measurement_t statusAt(const double& z) const //return also error in future
      {
	double x=position_.X()+angle_.X()*z;
	double y=position_.Y()+angle_.Y()*z;
	return Measurement_t(x,y);
      }
      
      void addMeasurement(TrackMeasurement& hit)
      {
	hits_.push_back(hit);
      }
      
      double chi2() 
      {
      };

      void fitTrack() 
      {
      };
      
      Measurement_t position_; //x,y
      MeasurementErrorMatrix_t positionError_;
      Measurement_t angle_; //alpha,beta angles in xz & yz
      MeasurementErrorMatrix_t angleError_;
      const TelescopeLayout& hodo_;
      std::vector<TrackMeasurement> hits_;
    };
   
    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:

    std::vector<Track> tracks_;
    //    TrackTree*     trackTree_;
};

DEFINE_PLUGIN(TrackReco);

#endif
