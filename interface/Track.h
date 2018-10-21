#include "Math/SMatrix.h"
#include "Math/SVector.h"

#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"
#include "Math/AxisAngle.h"
#include "Math/DisplacementVector3D.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "TObject.h"
#include "TNamed.h"

#include <iostream>

using namespace std;

using TrackParameters_t = ROOT::Math::SVector<double,4>; 
using TrackParametersCovarianceMatrix_t = ROOT::Math::SMatrix<double,4,4,ROOT::Math::MatRepSym<double,4> >;
using Measurement_t = ROOT::Math::SVector<double,2>;
using LocalCoord_t = ROOT::Math::SVector<double,2>;
using GlobalCoord_t = ROOT::Math::SVector<double,3>;
using MeasurementErrorMatrix_t = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >;
using RotationMatrix_t = ROOT::Math::SMatrix<double,3,3>;
using LocalRotationMatrix_t = ROOT::Math::SMatrix<double,2,2>;

namespace Tracking {
  class TrackLayer
  {
  public:
    TrackLayer()
      {
      };
  
  TrackLayer(const GlobalCoord_t& pos): 
    position_(pos),rotation_( ROOT::Math::SMatrixIdentity() )
      {
      };
    
  TrackLayer(const GlobalCoord_t& pos, const RotationMatrix_t& rot):
    position_(pos),rotation_(rot)
    {

    };

  TrackLayer(const GlobalCoord_t& pos, const double& zRot):
    position_(pos), rotation_(), zRotation_(zRot)
    {
      ROOT::Math::Rotation3D::Scalar zRotationAngle = zRotation_;
      ROOT::Math::RotationZ r_z(zRotationAngle);      
      ROOT::Math::Rotation3D rotation(r_z);
      std::vector<double> rot_components(9);
      rotation.GetComponents(rot_components.begin());
      rotation_.SetElements(rot_components.begin(),rot_components.end());
    };

    virtual ~TrackLayer()
      {
      };

    void setZRotation(const double& zRot)
    {
      zRotation_=zRot;
      ROOT::Math::Rotation3D::Scalar zRotationAngle = zRotation_;
      ROOT::Math::RotationZ r_z(zRotationAngle);      
      ROOT::Math::Rotation3D rotation(r_z);
      std::vector<double> rot_components(9);
      rotation.GetComponents(rot_components.begin());
      rotation_.SetElements(rot_components.begin(),rot_components.end());
    }

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

    inline bool measureX()
    {
      return measurementType_ & 0x1;
    }

    inline bool measureY()
    {
      return (measurementType_>>1) & 0x1 ;
    }

    /*
      inline void localCoordinates( GlobalCoord_t& loc, GlobalCoord_t& pos ) const;  
      { 
      } 
    */
      
    GlobalCoord_t position_; //axis origin in global Frame 
    RotationMatrix_t rotation_; //axis rotation matrix from localFrame to global Frame 
    double zRotation_;
    int measurementType_; //0=no measurement passive layer,1=x only,2=y only,3=x,y only

    ClassDef(TrackLayer, 3)
  };

  class TelescopeLayout : public TNamed
  {
  public:

  TelescopeLayout() :  TNamed()
      {
	layers_.clear();
      };

    TelescopeLayout(CfgManager& opts, string tagName);
      
    void addLayer(const TrackLayer& layer)
    {
      layers_.push_back(layer);
    }

    void SetName(const char *name)
    {
      fName = name;
    }

    inline void globalCoordinates( const int& k, const Measurement_t& loc, Measurement_t& pos ) const
    {
      layers_[k].globalCoordinates(loc, pos);
    }

    inline void globalCoordinates( const int& k, const MeasurementErrorMatrix_t& locErr, MeasurementErrorMatrix_t& posErr ) const
    {
      layers_[k].globalCoordinates(locErr, posErr);
    }

    void Print()
    {
      std::cout << "=== TELESCOPE GEOMETRY ===" << std::endl;
      int i=0;
      for (auto& layer : layers_)
	{
	  std::cout << "LAYER "<< i << " ===> [" << layer.position_ << "]:[" << layer.zRotation_ << "]" <<  std::endl;
	  ++i;
	}
    }

    std::vector<TrackLayer> layers_;

    ClassDef(TelescopeLayout, 4)
  };

  class TrackMeasurement : public TObject 
  {
  public:

  TrackMeasurement( Measurement_t position, MeasurementErrorMatrix_t positionErr,  const TelescopeLayout* hodo=NULL, int layer=0 ) : 
    localPosition_(position), localPositionError_ (positionErr), hodo_(hodo), layer_(layer)
    {
      calculateInverseVariance();
    }
      
  TrackMeasurement( double x, double y, const TelescopeLayout* hodo=NULL, int layer=0 ) :
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
      hodo_->globalCoordinates(layer_, localPosition_, pos);
    }

    inline void globalPositionError( MeasurementErrorMatrix_t& posErr ) const 
    {
      hodo_->globalCoordinates(layer_, localPositionError_, posErr);
    }

    inline void globalPositionErrorInverse( MeasurementErrorMatrix_t& posErr ) const 
    {
      hodo_->globalCoordinates(layer_, localPositionErrorInverse_, posErr);
    }

    Measurement_t localPosition_;
    MeasurementErrorMatrix_t localPositionError_; //covariance matrix
    MeasurementErrorMatrix_t localPositionErrorInverse_; //covariance matrix^-1
    const TelescopeLayout* hodo_;
    int layer_;
  };
     
  class LayerMeasurements: public TObject
  {
  public:
    LayerMeasurements() { hits_.clear(); };
    std::vector<TrackMeasurement> hits_;
  };
 
  class Track : public TObject
  {
  public:
  Track(TelescopeLayout* hodo=NULL) : trackPattern_(0), covarianceMatrixStatus_(0), hodo_(hodo)
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
      trackPattern_ |= 1 << hit.layer_;
    }

    void removeMeasurement(std::vector<TrackMeasurement>::iterator it)
    {
      int layer=it->layer_;
      hits_.erase(it);
      trackPattern_ &= 0 << layer;
    }

    void setTelescopeLayout(TelescopeLayout* hodo)
    {
      hodo_=hodo;
      for (auto& hit: hits_)
	hit.hodo_=hodo;
    }
    
    inline int nFreeParameters()
    {
      //count pixel hits as 2
      int nFreePar=0;
      for (auto& hit: hits_)
	{
	  if (hodo_->layers_[hit.layer_].measurementType_ == 3)
	    nFreePar+=2;
	  else if(hodo_->layers_[hit.layer_].measurementType_ > 0)
	    nFreePar+=1;
	}
      return nFreePar;
    }

      
    double chi2(const double* par=NULL); 

    double residual(const TrackMeasurement& hit,bool print=false);

    bool fitTrack();

    TrackParameters_t trackPar_; //x,y,alpha(xz angle),beta(yz angle)
    TrackParametersCovarianceMatrix_t trackParCov_; 
    unsigned int trackPattern_;
    int covarianceMatrixStatus_;
    TelescopeLayout* hodo_;
    std::vector<TrackMeasurement> hits_;
    bool fitAngle_;
  };

  class TrackContainer: public TObject
  {
  public:
    TrackContainer() { tracks_.clear(); };
    std::vector<Track> tracks_;
  };
  

}
