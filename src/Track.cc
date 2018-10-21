#include "interface/Track.h"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"
#include "Math/AxisAngle.h"
#include "Math/DisplacementVector3D.h"


double Tracking::Track::residual(const TrackMeasurement& hit,bool print)
{
  Measurement_t pos;
  MeasurementErrorMatrix_t posErrInverse;
  hit.globalPosition(pos);
  hit.globalPositionErrorInverse(posErrInverse);
  Measurement_t delta= pos-this->statusAt((this->hodo_->layers_[hit.layer_].position_(2))); //calculation of the residual
  if (print)
    {
      std::cout << "HIT POS " << pos << std::endl;
      std::cout << "TRACK STATUS " << this->statusAt((this->hodo_->layers_[hit.layer_].position_(2))) << std::endl;
      std::cout << "RESIDUAL " << ROOT::Math::Dot(delta,posErrInverse*delta) << std::endl;
    }
  return ROOT::Math::Dot(delta,posErrInverse*delta); // delta^T * (V^-1) * delta
}

double Tracking::Track::chi2(const double* par)
{
  if (par) //--- set the actual fit parameters
    {
      trackPar_(0)=par[0];
      trackPar_(1)=par[1];
      if (fitAngle_)
	{
	  trackPar_(2)=par[2];
	  trackPar_(3)=par[3];
	}
    }
  
  double chi2=0;
  for (auto& hit: Tracking::Track::hits_)
      chi2+=residual(hit);

  return chi2;
}


bool Tracking::Track::fitTrack()
{
  //---setup minimization
  int nFitParameters = fitAngle_? 4:2;
  ROOT::Math::Functor chi2(this, &Tracking::Track::chi2, nFitParameters);
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(1000000);
  minimizer->SetTolerance(1e-6);
  minimizer->SetPrintLevel(0);
  minimizer->SetFunction(chi2);
  minimizer->SetLimitedVariable(0, "X", trackPar_(0), 1E-6, -30, 30);
  minimizer->SetLimitedVariable(1, "Y", trackPar_(1), 1E-6, -30, 30);
  if (fitAngle_)
    {
      minimizer->SetLimitedVariable(2, "alpha",trackPar_(2), 1E-6,  -0.1, 0.1);
      minimizer->SetLimitedVariable(3, "beta" ,trackPar_(3), 1E-6,  -0.1, 0.1);
    }

  //---fit
  minimizer->Minimize();

  //---fill best fit par
  for (int i=0;i<nFitParameters;++i)
    trackPar_(i) = minimizer->X()[i];

  //---get covariance matrix
  covarianceMatrixStatus_ = minimizer->CovMatrixStatus();
  double covariances[nFitParameters * nFitParameters];
  minimizer->GetCovMatrix(covariances);

  //---fill covariance matrix
  for (int i=0;i<nFitParameters * nFitParameters;++i)
    if ( i/nFitParameters >= i%nFitParameters)
      trackParCov_(i/nFitParameters,i%nFitParameters)=covariances[i];
  
  delete minimizer;        

  return true;
}

Tracking::TelescopeLayout::TelescopeLayout(CfgManager& opts,string tagName)
{
  //---inputs---
  std::vector<string> layers = opts.GetOpt<vector<string> >(tagName+".layers");
  
  for (auto& layer: layers)
    {
      int measurementType=opts.GetOpt<int>(layer+".measurementType");
      std::vector<double> position=opts.GetOpt<vector<double> >(layer+".position");
      if (position.size() != 3)
	std::cout << "ERROR: Expecting a vector of size 3 for the layer position" << std::endl;
      GlobalCoord_t layerPos;
      layerPos.SetElements(position.begin(),position.end());
      
      ROOT::Math::Rotation3D::Scalar zRotationAngle = opts.GetOpt<ROOT::Math::Rotation3D::Scalar>(layer+".zRotationAngle");

      Tracking::TrackLayer aLayer(layerPos,zRotationAngle);
      aLayer.measurementType_=measurementType;
      addLayer(aLayer);      
    }
}
