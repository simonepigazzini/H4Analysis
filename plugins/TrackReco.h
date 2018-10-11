#ifndef __TRACK_RECO__
#define __TRACK_RECO__

#include "interface/PluginBase.h"
#include "interface/TrackTree.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"

using namespace std;

using TrackParameters_t = ROOT::Math::SVector<double,4>; 
using TrackParametersCovarianceMatrix_t = ROOT::Math::SMatrix<double,4,4,ROOT::Math::MatRepSym<double,4> >;
using Measurement_t = ROOT::Math::SVector<double,2>;
using LocalCoord_t = ROOT::Math::SVector<double,2>;
using GlobalCoord_t = ROOT::Math::SVector<double,3>;
using MeasurementErrorMatrix_t = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >;
using RotationMatrix_t = ROOT::Math::SMatrix<double,3,3>;
using LocalRotationMatrix_t = ROOT::Math::SMatrix<double,2,2>;

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
      position_(pos),rotation_( ROOT::Math::SMatrixIdentity() )
      {
      };

    TrackLayer(const GlobalCoord_t& pos, const RotationMatrix_t& rot):
      position_(pos),rotation_(rot)
      {
      };
      
      inline void globalCoordinates( const Measurement_t& loc, Measurement_t& pos ) const
      {
	pos= rotation_.Sub<LocalRotationMatrix_t>(0,0)*loc + position_.Sub<LocalCoord_t>(0);
      };

      inline void globalCoordinates( const MeasurementErrorMatrix_t& err, MeasurementErrorMatrix_t& globalError ) const
      {
	//LocalRotationMatrix_t r=rotation_.Sub<LocalRotationMatrix_t>(0,0);
	//globalError=rotation_.Sub<LocalRotationMatrix_t>(0,0)*err*ROOT::Math::Transpose(rotation_.Sub<LocalRotationMatrix_t>(0,0));
	globalError= err;
      };
      
      /*
	inline void localCoordinates( GlobalCoord_t& loc, GlobalCoord_t& pos ) const;  
	 { 
	 } 
      */
      
      GlobalCoord_t position_; //axis origin in global Frame 
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

      inline void globalCoordinates( const int& k, const Measurement_t& loc, Measurement_t& pos ) const
      {
	layers_[k].globalCoordinates(loc, pos);
      }

      inline void globalCoordinates( const int& k, const MeasurementErrorMatrix_t& locErr, MeasurementErrorMatrix_t& posErr ) const
      {
	layers_[k].globalCoordinates(locErr, posErr);
      }
            
      std::vector<TrackLayer> layers_;
    };
    
    class TrackMeasurement {
    public:
    TrackMeasurement( double x, double y, const TelescopeLayout& hodo, int layer ) :
      localPosition_(x,y), localPositionError_(), hodo_(hodo), layer_(layer)
      {
	//check that layer is within the size of telescope...
      }
      
      ~TrackMeasurement() {};
      
      inline void setVarianceX(const double& sigma2x) 
      {
	localPositionError_(0,0)=sigma2x;
      }

      inline void setVarianceY(const double& sigma2y) 
      {
	localPositionError_(1,1)=sigma2y;
      }

      inline void calculateInverseVariance()
      {
	int ifail;
	localPositionErrorInverse_ = localPositionError_.InverseFast(ifail);
	if (ifail)
	  cout << "[TrackMeasurement]::[ERROR]::Matrix Inversion failed" << endl;
      }

      inline void globalPosition( Measurement_t& pos ) const 
      {
	hodo_.globalCoordinates(layer_, localPosition_, pos);
      }

      inline void globalPositionError( MeasurementErrorMatrix_t& posErr ) const 
      {
	hodo_.globalCoordinates(layer_, localPositionError_, posErr);
      }

      inline void globalPositionErrorInverse( MeasurementErrorMatrix_t& posErr ) const 
      {
	hodo_.globalCoordinates(layer_, localPositionErrorInverse_, posErr);
      }

      Measurement_t localPosition_;
      MeasurementErrorMatrix_t localPositionError_; //covariance matrix
      MeasurementErrorMatrix_t localPositionErrorInverse_; //covariance matrix^-1
      const TelescopeLayout& hodo_;
      int layer_;
    };
  
    class Track {
    public:
    Track(const TelescopeLayout& hodo) : covarianceMatrixStatus_(0), hodo_(hodo)
      {
	hits_.clear();
      };
      
      ~Track() {};
      
      inline Measurement_t statusAt(const double& z) const //return also error in future
      {
	return trackPar_.Sub<Measurement_t>(0) + trackPar_.Sub<Measurement_t>(2) * z;
      }
      
      void addMeasurement(TrackMeasurement& hit)
      {
	hits_.push_back(hit);
      }
      
      double chi2(const double* par=NULL); 

      bool fitTrack();

      TrackParameters_t trackPar_; //x,y,alpha(xz angle),beta(yz angle)
      TrackParametersCovarianceMatrix_t trackParCov_; 
      int covarianceMatrixStatus_;
      const TelescopeLayout& hodo_;
      std::vector<TrackMeasurement> hits_;
      bool fitAngle_;
    };
   
    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:

    std::vector<Track> tracks_;
    TelescopeLayout hodo_;
    TrackTree*     trackTree_;
};

DEFINE_PLUGIN(TrackReco);

#endif
