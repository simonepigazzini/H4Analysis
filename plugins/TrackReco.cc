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
  for (auto& hit: TrackReco::Track::hits_)
    {
      //calculate the residual
      Measurement_t pos;
      MeasurementErrorMatrix_t posErrInverse;
      hit.globalPosition(pos);
      hit.globalPositionErrorInverse(posErrInverse);
      Measurement_t delta= pos-this->statusAt((this->hodo_.layers_[hit.layer_].position_(2))); //calculation of the residual
      chi2+=ROOT::Math::Dot(delta,posErrInverse*delta); // delta^T * (V^-1) * delta
    }
  return chi2;
}

bool TrackReco::Track::fitTrack()
{
  //---setup minimization
  int nFitParameters = fitAngle_? 4:2;
  ROOT::Math::Functor chi2(this, &TrackReco::Track::chi2, nFitParameters);
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

//----------Begin-------------------------------------------------------------------------
bool TrackReco::Begin(CfgManager& opts, uint64* index)
{
  //---create a position tree
  bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
    opts.GetOpt<bool>(instanceName_+".storeTree") : true;
  
  //---create a position tree
  string treeName = opts.OptExist(instanceName_+".treeName") ?
    opts.GetOpt<string>(instanceName_+".treeName") : "track_tree";

  RegisterSharedData(new TTree(treeName.c_str(), treeName.c_str()), treeName.c_str(), storeTree);
  
  trackTree_ = new TrackTree(index, (TTree*)data_.back().obj);
  trackTree_->Init();
  
  //  Now build the hodoscope - just an example
  GlobalCoord_t layer0Pos(0.,0.,0.);
  GlobalCoord_t layer1Pos(0.,0.,100);
  
  GlobalCoord_t layer2Pos(0.,0.,500.);
  RotationMatrix_t layer2Rot;
  layer2Rot=ROOT::Math::SMatrixIdentity();
  // double angle=-0.05;
  // layer2Rot(0,0)=cos(angle);
  // layer2Rot(0,1)=sin(angle);
  // layer2Rot(1,0)=-sin(angle);
  // layer2Rot(1,1)=cos(angle);

  GlobalCoord_t layer3Pos(0.,0.,1000);
  GlobalCoord_t layer4Pos(0.,0.,1100);
  
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
  trackTree_->Clear();

  TRandom3 r(0);
  //Get track measurements
  for (int i=0;i<1000;++i)
    {
      TrackMeasurement aHit_0(r.Gaus(0.3,0.2),0.,hodo_,0); //X measurement
      aHit_0.setVarianceX(0.2*0.2);
      aHit_0.setVarianceY(9999.);
      aHit_0.calculateInverseVariance(); //important to be done after filling/setting the variance
      TrackMeasurement aHit_1(0.,r.Gaus(0.5,0.2),hodo_,1); //Y measurement
      aHit_1.setVarianceX(9999.);
      aHit_1.setVarianceY(0.2*0.2);
      aHit_1.calculateInverseVariance();
      TrackMeasurement aHit_2(r.Gaus(1.3,0.02),r.Gaus(0.5,0.02),hodo_,2); //X,Y measurement
      aHit_2.setVarianceX(0.02*0.02);
      aHit_2.setVarianceY(0.02*0.02);
      aHit_2.calculateInverseVariance();
      TrackMeasurement aHit_3(r.Gaus(2.3,0.2),0.,hodo_,3); //X measurement
      aHit_3.setVarianceX(0.2*0.2);
      aHit_3.setVarianceY(9999.);
      aHit_3.calculateInverseVariance();
      TrackMeasurement aHit_4(0.,r.Gaus(0.5,0.2),hodo_,4); //Y measurement
      aHit_4.setVarianceX(9999.);
      aHit_4.setVarianceY(0.2*0.2);
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

      tracks_.push_back(aTrack);
    }
  //---fill output tree
  trackTree_->n_tracks=tracks_.size();
  for (auto& aTrack: tracks_)
    {
      TrackPar par;
      par.value.assign(aTrack.trackPar_.Array(),aTrack.trackPar_.Array()+4);
      par.covariance.assign(aTrack.trackParCov_.Array(),aTrack.trackParCov_.Array()+10);
      trackTree_->fitResult.push_back(par);
      trackTree_->fitStatus.push_back(aTrack.covarianceMatrixStatus_);
      trackTree_->trackHits.push_back(aTrack.hits_.size());
      trackTree_->trackChi2.push_back(aTrack.chi2());
    }
  trackTree_->Fill();
  return true;
}
