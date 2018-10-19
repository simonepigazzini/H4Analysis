#include "AlignHodo.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"
#include "Math/AxisAngle.h"
#include "Math/DisplacementVector3D.h"

//----------Begin-------------------------------------------------------------------------
bool AlignHodo::Begin(CfgManager& opts, uint64* index)
{  
  //---inputs---
  srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstance");
  
  //  hodo_=0;
  tracks_.tracks_.clear();
  return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool AlignHodo::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
  //---load tracks from trackReco
  string objectName="tracks";
  auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_tracks", "", false);

  Tracking::TrackContainer* tracks;

  if(shared_data.size() != 0)
    tracks = (Tracking::TrackContainer*)shared_data.at(0).obj;
  else
    cout << "[AlignHodo::" << instanceName_ << "]: " << srcInstance_+"_tracks not found" << endl; 

  //---select good tracks
  for (auto& track : tracks->tracks_)
    {
      if (hodo_.layers_.size() == 0)
	hodo_=*track.hodo_;

      if (track.hits_.size()>3) 
	tracks_.tracks_.push_back(track);
    }
  
  return true;
}

double AlignHodo::globalChi2(const double* par)
{
  double chi2=0;

  Tracking::TelescopeLayout alignedHodo(hodo_); //reset hodoscope
  std::vector<Tracking::TrackLayer> originalLayers;
  std::vector<Tracking::TrackLayer> alignedLayers;
  originalLayers.assign(alignedHodo.layers_.begin(),alignedHodo.layers_.end());
  alignedLayers.assign(alignedHodo.layers_.begin(),alignedHodo.layers_.end());

  if (par) //--- set the actual fit parameters
    {
      for (int i=0;i<(hodo_.layers_.size()-2);++i) //keep the first 2 layers as reference
	{
	  //--- set position and rotations
	  for (int icoord=0;icoord<3;++icoord) 
	    alignedLayers[i+2].position_(icoord)=par[i*4+icoord];
	  
	  ROOT::Math::Rotation3D::Scalar zRotationAngle = par[i*4+3];
	  ROOT::Math::RotationZ r_z(zRotationAngle);      
	  ROOT::Math::Rotation3D rotation(r_z);
	  std::vector<double> rot_components(9);
	  rotation.GetComponents(rot_components.begin());
	  alignedLayers[i+2].rotation_.SetElements(rot_components.begin(),rot_components.end());
	}
    }

  for (auto& track: tracks_.tracks_)
    {
      int iHit=0;
      for (auto it=track.hits_.begin(); it != track.hits_.end(); ++it)
      	{
      	  if ( it->layer_ > 1)
      	    {
      	      Tracking::TrackMeasurement aHit(*it);
      	      Tracking::Track aTrack=track;
	      alignedHodo.layers_[it->layer_]=alignedLayers[it->layer_]; //put misaligned layer (do this layer by layer)
	      aHit.hodo_=&alignedHodo;
	      aTrack.setTelescopeLayout(&alignedHodo);
	      aTrack.fitAngle_=false;
	      aTrack.removeMeasurement(aTrack.hits_.begin()+iHit);
	      aTrack.fitTrack(); //now refit without the hit
      	      chi2+=aTrack.residual(aHit); //take output chi2
	      alignedHodo.layers_[it->layer_]=originalLayers[it->layer_]; //restore layer original position  
      	    }
	  ++iHit;
    	}
    }

  std::cout << "***" << chi2/tracks_.tracks_.size() << std::endl;
  return chi2;
}

void AlignHodo::minimize()
{
  //---setup minimization
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(1000000);
  minimizer->SetTolerance(1e-6);
  minimizer->SetPrintLevel(1);
  //---setup variables and ranges
  int nFitParameters = (hodo_.layers_.size()-2)*4;
  for (int i=0;i<(hodo_.layers_.size()-2);++i) //keep the first 2 layers as reference
    {
      GlobalCoord_t layerPos=hodo_.layers_[i+2].position_;
      minimizer->SetLimitedVariable(i*4, Form("X_%d",i+2), layerPos(0), 1E-6, -30, 30);
      minimizer->SetLimitedVariable(i*4+1, Form("Y_%d",i+2), layerPos(1), 1E-6, -30, 30);
      minimizer->SetLimitedVariable(i*4+2, Form("Z_%d",i+2), layerPos(2), 1E-6, 0, 3000);
      if (hodo_.layers_[i+2].measurementType_==3)
      	minimizer->SetLimitedVariable(i*4+3, Form("zRot_%d",i+2), 0, 1E-6, -30, 30);
      else
      	minimizer->SetFixedVariable(i*4+3, Form("zRot_%d",i+2), 0);
    }

  std::cout << "START MINIMIZATION" << std::endl;
  //---fit
  ROOT::Math::Functor chi2(this, &AlignHodo::globalChi2, nFitParameters);
  minimizer->SetFunction(chi2);
  minimizer->Minimize();

  // //---fill best fit par

  for (int i=0;i<(hodo_.layers_.size()-2);++i) //keep the first 2 layers as reference
    {
      //--- set position and rotations
      for (int icoord=0;icoord<3;++icoord) 
	hodo_.layers_[i+2].position_(icoord)=minimizer->X()[i*4+icoord];
      
      ROOT::Math::Rotation3D::Scalar zRotationAngle = minimizer->X()[i*4+3];
      ROOT::Math::RotationZ r_z(zRotationAngle);      
      ROOT::Math::Rotation3D rotation(r_z);
      std::vector<double> rot_components(9);
      rotation.GetComponents(rot_components.begin());
      hodo_.layers_[i+2].rotation_.SetElements(rot_components.begin(),rot_components.end());
    }

  // //---get covariance matrix
  // covarianceMatrixStatus_ = minimizer->CovMatrixStatus();
  // double covariances[nFitParameters * nFitParameters];
  // minimizer->GetCovMatrix(covariances);

  // //---fill covariance matrix
  // for (int i=0;i<nFitParameters * nFitParameters;++i)
  //   if ( i/nFitParameters >= i%nFitParameters)
  //     trackParCov_(i/nFitParameters,i%nFitParameters)=covariances[i];
  
  delete minimizer;        
}

bool AlignHodo::End(CfgManager& opts)
{
  minimize();
  TFile *f=new TFile("align.root","RECREATE");
  TTree* t_ali=new TTree("alignHodo","alignHodo");
  t_ali->Branch("alignedHodo",&hodo_);
  t_ali->Fill();
  t_ali->Write();
  f->Close();
  return true;
}
