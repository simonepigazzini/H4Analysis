#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TObject.h"

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
        TelescopeLayer (const GlobalCoord_t& pos):
            position_(pos),rotation_( ROOT::Math::SMatrixIdentity() )
            {
            };
    
        TelescopeLayer (const GlobalCoord_t& pos, const RotationMatrix_t& rot):
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

        inline bool measureX()
            {
                return measurementType_.find("X") != string::npos;
            }

        inline bool measureY()
            {
                return measurementType_.find("Y") != string::npos;
            }

        /*
          inline void localCoordinates( GlobalCoord_t& loc, GlobalCoord_t& pos ) const;  
          { 
          } 
        */
      
        GlobalCoord_t    position_; //axis origin in global Frame 
        RotationMatrix_t rotation_; //axis rotation matrix from localFrame to global Frame 
        string           measurementType_;
    };

    //---Telescope layout: layer container (z relative position are recorded here
    class TelescopeLayout
    {
    public:
        TelescopeLayout()
            {
                layers_.clear();
            };
      
        void addLayer(const TelescopeLayer & layer)
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
                    std::cout << "LAYER "<< i << " ===> [" << layer.position_ << "]" << std::endl;
                    ++i;
                }
            }

        std::vector<TelescopeLayer> layers_;
    };

    //---Single layer, single hit measurement
    class TrackHit : public TObject 
    {
    public:

        TrackHit( Measurement_t position, MeasurementErrorMatrix_t positionErr, const TelescopeLayout* tLayout=NULL, int layer=0 ) : 
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

        Measurement_t            localPosition_;
        MeasurementErrorMatrix_t localPositionError_; //covariance matrix
        MeasurementErrorMatrix_t localPositionErrorInverse_; //covariance matrix^-1
        const TelescopeLayout*   tLayout_;
        int                      layer_;
    };

    //---Single layer measurement: this is a collection of all hits in a given layer
    class LayerHits: public TObject
    {
    public:
        LayerHits() { hits_.clear(); };
        std::vector<TrackHit> hits_;
    };

    //---Track object (collection of multiple layer measurements):
    //   - perform track fit
    class Track : public TObject
    {
    public:
        Track(const TelescopeLayout* tLayout=NULL) : trackPattern_(0), covarianceMatrixStatus_(0), tLayout_(tLayout)
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
      
        double chi2(const double* par=NULL); 

        double residual(const TrackHit& hit,bool print=false);

        bool fitTrack();

        TrackParameters_t                 trackPar_; //x,y,alpha(xz angle),beta(yz angle)
        TrackParametersCovarianceMatrix_t trackParCov_; 
        unsigned int                      trackPattern_;
        int                               covarianceMatrixStatus_;
        const TelescopeLayout*            tLayout_;
        std::vector<TrackHit>             hits_;
        bool                              fitAngle_;
    };
}
