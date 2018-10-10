#include "TrackReco.h"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"
#include "Math/AxisAngle.h"
#include "Math/DisplacementVector3D.h"
#include "TRandom3.h"

double TrackReco::Track::chi2(const double* par)
{
  if (par) //set the actual fit parameters
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
   int i=0;
   for (auto& hit: TrackReco::Track::hits_)
     {
       //calculate the residual
       Measurement_t pos;
       MeasurementErrorMatrix_t posErrInverse;
       hit.globalPosition(pos);
       hit.globalPositionErrorInverse(posErrInverse);
       Measurement_t delta= pos-this->statusAt((this->hodo_.layers_[hit.layer_].position_(2))); //calculation of the residual
       chi2+=ROOT::Math::Dot(delta,posErrInverse*delta); // delta^T * (V^-1) * delta
       ++i;
     }
   return chi2;
}

bool TrackReco::Track::fitTrack()
{
  //---setup minimization
  int nFreeParameters = fitAngle_? 4:2;
  ROOT::Math::Functor chi2(this, &TrackReco::Track::chi2, nFreeParameters);
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

  for (int i=0;i<nFreeParameters;++i)
    trackPar_(i) = minimizer->X()[i];

  covarianceMatrixStatus_ = minimizer->CovMatrixStatus();
  double covariances[nFreeParameters * nFreeParameters];
  minimizer->GetCovMatrix(covariances);

  for (int i=0;i<nFreeParameters * nFreeParameters;++i)
    if ( i/nFreeParameters >= i%nFreeParameters)
      trackParCov_(i/nFreeParameters,i%nFreeParameters)=covariances[i];
  
  delete minimizer;        

  return true;
}

//----------Begin-------------------------------------------------------------------------
bool TrackReco::Begin(CfgManager& opts, uint64* index)
{
    //---create a position tree
//    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
//        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

//  Now build the hodoscope - just an example for now
    GlobalCoord_t layer0Pos(0.,0.,0.);
    GlobalCoord_t layer1Pos(0.,0.,100);
    GlobalCoord_t layer2Pos(0.,0.,0.);
    RotationMatrix_t layer2Rot;
    layer2Rot=ROOT::Math::SMatrixIdentity();
    // double angle=-0.05;
    // layer2Rot(0,0)=cos(angle);
    // layer2Rot(0,1)=sin(angle);
    // layer2Rot(1,0)=-sin(angle);
    // layer2Rot(1,1)=cos(angle);
    
    GlobalCoord_t layer3Pos(0.,0.,2000);
    GlobalCoord_t layer4Pos(0.,0.,2100);

    TrackLayer layer_0(layer0Pos);
    TrackLayer layer_1(layer1Pos);
    TrackLayer layer_2(layer2Pos,layer2Rot);
    TrackLayer layer_3(layer3Pos);
    TrackLayer layer_4(layer4Pos);

    hodo_.addLayer(layer_0);
    hodo_.addLayer(layer_1);
    hodo_.addLayer(layer_2);
    hodo_.addLayer(layer_3);
    hodo_.addLayer(layer_4);

    std::cout << "BUILT HODO" << std::endl;
    int i=0;
    for (auto& layer : hodo_.layers_)
      {
	cout << "LAYER "<< i << ":" << layer.position_ << std::endl;
	++i;
      }
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool TrackReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{

    tracks_.clear();

    TRandom3 r(0);
    //Get track measurements
    for (int i=0;i<100;++i)
      {
	std::cout << "=======================" << std::endl;
	TrackMeasurement aHit_0(r.Gaus(0.3,0.2),0.,hodo_,0);
	aHit_0.setVarianceX(0.04);
	aHit_0.setVarianceY(9999.);
	aHit_0.calculateInverseVariance();
	TrackMeasurement aHit_1(0.,r.Gaus(0.6,0.2),hodo_,1);
	aHit_1.setVarianceX(9999.);
	aHit_1.setVarianceY(0.04);
	aHit_1.calculateInverseVariance();
	TrackMeasurement aHit_2(r.Gaus(0.3,0.02),r.Gaus(0.6,0.02),hodo_,2);
	aHit_2.setVarianceX(0.0004);
	aHit_2.setVarianceY(0.0004);
	aHit_2.calculateInverseVariance();
	TrackMeasurement aHit_3(r.Gaus(0.3,0.2),0.,hodo_,3);
	aHit_3.setVarianceX(0.04);
	aHit_3.setVarianceY(9999.);
	aHit_3.calculateInverseVariance();
	TrackMeasurement aHit_4(0.,r.Gaus(0.6,0.2),hodo_,4);
	aHit_4.setVarianceX(9999.);
	aHit_4.setVarianceY(0.04);
	aHit_4.calculateInverseVariance();
	
	Track aTrack(hodo_);
	aTrack.addMeasurement(aHit_0);
	aTrack.addMeasurement(aHit_1);
	aTrack.addMeasurement(aHit_2);
	aTrack.addMeasurement(aHit_3);
	aTrack.addMeasurement(aHit_4);
	
	//fit track
	aTrack.fitAngle_ = true;
	aTrack.fitTrack();

	std::cout << "BEST FIT " << aTrack.trackPar_(0) << "," << aTrack.trackPar_(1) ;
	if (aTrack.fitAngle_)
	  std::cout << ":" << aTrack.trackPar_(2) << "," << aTrack.trackPar_(3) ;
	std::cout << std::endl;

	std::cout << "CHI2 " << aTrack.chi2() << std::endl;
	std::cout << "ERR X " << sqrt(aTrack.trackParCov_(0,0)) << std::endl;
	std::cout << "ERR Y " << sqrt(aTrack.trackParCov_(1,1)) << std::endl;
	if (aTrack.fitAngle_)
	  {
	    std::cout << "ERR ALPHA " << sqrt(aTrack.trackParCov_(2,2)) << std::endl;
	    std::cout << "ERR BETA " << sqrt(aTrack.trackParCov_(3,3)) << std::endl;
	  }
      }
    //---fill output tree
    return true;
}
