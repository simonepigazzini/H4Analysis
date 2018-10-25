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
    class TelescopeLayer
    {
    public:
        TelescopeLayer()
            {
            };
  
        TelescopeLayer(const GlobalCoord_t& pos): 
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
        string measurementType_; // "", "X", "Y", "XY"

        ClassDef(TelescopeLayer, 3)
    };

    //---Telescope layout: contain the structure of the telescope and reference to its layers
    class TelescopeLayout : public TNamed
    {
    public:

        TelescopeLayout() :  TNamed()
            {
                layers_.clear();
            };

        TelescopeLayout(CfgManager& opts, string tagName);
      
        void addLayer(const TelescopeLayer& layer)
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

        std::vector<TelescopeLayer> layers_;

        ClassDef(TelescopeLayout, 5)
    };

    //---Single track hit in a telescope layer
    class TrackHit : public TObject 
    {
    public:

        TrackHit( Measurement_t position, MeasurementErrorMatrix_t positionErr,  const TelescopeLayout* tLayout=NULL, int layer=0 ) : 
            localPosition_(position), localPositionError_ (positionErr), tLayout_(tLayout), layer_(layer)
            {
                calculateInverseVariance();
            }
      
        TrackHit( double x, double y, const TelescopeLayout* tLayout=NULL, int layer=0 ) :
            localPosition_(x,y), localPositionError_(), tLayout_(tLayout), layer_(layer)
            {
                //check that layer is within the size of telescope...
            }
      
        ~TrackHit() {};
      
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
                    cout << "[TrackHit]::[ERROR]::Matrix Inversion failed" << endl;
            }

        inline void globalPosition( Measurement_t& pos ) const 
            {
                tLayout_->globalCoordinates(layer_, localPosition_, pos);
            }

        inline void globalPositionError( MeasurementErrorMatrix_t& posErr ) const 
            {
                tLayout_->globalCoordinates(layer_, localPositionError_, posErr);
            }

        inline void globalPositionErrorInverse( MeasurementErrorMatrix_t& posErr ) const 
            {
                tLayout_->globalCoordinates(layer_, localPositionErrorInverse_, posErr);
            }

        Measurement_t localPosition_;
        MeasurementErrorMatrix_t localPositionError_; //covariance matrix
        MeasurementErrorMatrix_t localPositionErrorInverse_; //covariance matrix^-1
        const TelescopeLayout* tLayout_;
        int layer_;
    };

    //---Collection of hits (helper class)
    class LayerHits: public TObject
    {
    public:
        LayerHits() { hits_.clear(); };
        std::vector<TrackHit> hits_;
    };

    //---Track object:
    //   - hits container
    //   - track fitting
    class Track : public TObject
    {
    public:
        Track(TelescopeLayout* tLayout=NULL) : trackPattern_(0), covarianceMatrixStatus_(0), tLayout_(tLayout)
            {
                hits_.clear();
            };
      
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

        TrackParameters_t trackPar_; //x,y,alpha(xz angle),beta(yz angle)
        TrackParametersCovarianceMatrix_t trackParCov_; 
        unsigned int trackPattern_;
        int covarianceMatrixStatus_;
        TelescopeLayout* tLayout_;
        std::vector<TrackHit> hits_;
        bool fitAngle_;
    };

    //---simple collection of tracks
    class TrackContainer: public TObject
    {
    public:
        TrackContainer() { tracks_.clear(); };
        std::vector<Track> tracks_;
    };
  

}
