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
    //---Single telescope layer
    class TelescopeLayer
    {
    public:
        TelescopeLayer()
            {
            };
  
        telescopeLayer(const GlobalCoord_t& pos): 
            position_(pos),rotation_( ROOT::Math::SMatrixIdentity() )
            {
            };
    
        TelescopeLayer(const GlobalCoord_t& pos, const RotationMatrix_t& rot):
            position_(pos),rotation_(rot)
            {

            };

        TelescopeLayer(const GlobalCoord_t& pos, const double& zRot):
            position_(pos), rotation_(), zRotation_(zRot)
            {
                ROOT::Math::Rotation3D::Scalar zRotationAngle = zRotation_;
                ROOT::Math::RotationZ r_z(zRotationAngle);      
                ROOT::Math::Rotation3D rotation(r_z);
                std::vector<double> rot_components(9);
                rotation.GetComponents(rot_components.begin());
                rotation_.SetElements(rot_components.begin(),rot_components.end());
            };

        virtual ~TelescopeLayer()
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
                return measurementType_.find("X") != string::npos;
            }

        inline bool measureY()
            {
                return measurementType_.find("Y") != string::npos;
            }
      
        GlobalCoord_t position_; //axis origin in global Frame 
        RotationMatrix_t rotation_; //axis rotation matrix from localFrame to global Frame 
        double zRotation_;
        string measurementType_; //0=no measurement passive layer,1=x only,2=y only,3=x,y only

        ClassDef(TelescopeLayer, 3)
    };

    //---Telescope layout record
    class TelescopeLayout : public TObject
    {
    public:

        TelescopeLayout()
            {
                layers_.clear();
            };

        TelescopeLayout(CfgManager& opts, string tagName);
      
        void addLayer(const TelescopeLayer& layer)
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


        std::vector<TelescopeLayer> layers_;

        ClassDef(TelescopeLayout, 3)
    };

    //---Single track hit in one telescope layer
    class TrackHit : public TObject 
    {
    public:

        TrackHit( Measurement_t position, MeasurementErrorMatrix_t positionErr,  const TelescopeLayout* tLayout=NULL, int layer=0 ) : 
            localPosition_(position), localPositionError_ (positionErr), tLayout_(tLayout), layer_(layer)
            {
                calculateInverseVariance();
            }
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

    //---All the hits reconded in one layer
    class LayerHits: public TObject
    {
    public:
        LayerHits() { hits_.clear(); };
        std::vector<TrackHit> hits_;
    };

    //---Track class:
    //   - hits container
    //   - track fitting
    class Track : public TObject
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
      
        ~Track() {};
      
        inline Measurement_t statusAt(const double& z) const //return also error in future
            {
                return trackPar_.Sub<Measurement_t>(0) + trackPar_.Sub<Measurement_t>(2) * z;
            }
      
        void addMeasurement(TrackHit& hit)
            {
                hits_.push_back(hit);
                trackPattern_ |= 1 << hit.layer_;
            }

        void removeMeasurement(std::vector<TrackHit>::iterator it)
            {
                int layer=it->layer_;
                hits_.erase(it);
                trackPattern_ &= 0 << layer;
            }

        void setTelescopeLayout(TelescopeLayout* tLayout)
            {
                tLayout_=tLayout;
                for (auto& hit: hits_)
                    hit.tLayout_=tLayout;
            }

        inline int nFreeParameters()
            {
                //count pixel hits as 2
                int nFreePar=0;
                for (auto& hit: hits_)
                {
                    if (tLayout_->layers_[hit.layer_].measurementType_ == "XY")
                        nFreePar+=2;
                    else if(tLayout_->layers_[hit.layer_].measurementType_.size() > 0)
                        nFreePar+=1;
                }
                return nFreePar;
            }

      
        double chi2(const double* par=NULL); 

        double residual(const TrackHit& hit,bool print=false);

        bool fitTrack();

        TrackParameters_t                 trackPar_; //x,y,alpha(xz angle),beta(yz angle)
        TrackParametersCovarianceMatrix_t trackParCov_; 
        unsigned int                      trackPattern_;
        int                               covarianceMatrixStatus_;
        TelescopeLayout*                  tLayout_;
        std::vector<TrackHit>             hits_;
        bool                              fitAngle_;
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
    public:
        TrackContainer() { tracks_.clear(); };
        std::vector<Track> tracks_;
    };
}
