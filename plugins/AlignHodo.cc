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

#include <time.h>

//----------Begin-------------------------------------------------------------------------
bool AlignHodo::Begin(CfgManager& opts, uint64* index)
{  
  //---inputs---
  srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstance");
  alignZ_ = opts.OptExist(instanceName_+".alignZ") ? opts.GetOpt<bool>(instanceName_+".alignZ") : true;

  hodo_=NULL;
  tracks_.tracks_.clear();
  return true;
}

//----------Begin-------------------------------------------------------------------------
bool AlignHodo::BeginLoop(int iLoop, CfgManager& opts)
{  
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
      if (!hodo_)
	hodo_=track.hodo_;

      if (track.hits_.size()>3) 
	tracks_.tracks_.push_back(track);
    }
  
  return true;
}

double AlignHodo::globalChi2(const double* par)
{
  double chi2=0;

  Tracking::TelescopeLayout alignedHodo(*hodo_); //reset hodoscope
  std::vector<Tracking::TrackLayer> originalLayers;
  std::vector<Tracking::TrackLayer> alignedLayers;
  originalLayers.assign(alignedHodo.layers_.begin(),alignedHodo.layers_.end());
  alignedLayers.assign(alignedHodo.layers_.begin(),alignedHodo.layers_.end());

  if (par) //--- set the actual fit parameters
    {
      for (int i=0;i<(hodo_->layers_.size()-2);++i) //keep the first 2 layers as reference
	{
	  //--- set position and rotations
	  for (int icoord=0;icoord<3;++icoord) 
	    alignedLayers[i+2].position_(icoord)=par[i*4+icoord];
	  alignedLayers[i+2].setZRotation(par[i*4+3]);

	  alignedHodo.layers_[i+2]=alignedLayers[i+2]; //put misaligned layer in layout 
	}
    }

  for (auto& track: tracks_.tracks_)
    {
      int iHit=0;
      for (auto it=track.hits_.begin(); it != track.hits_.end(); ++it)
      	{
      	  if ( it->layer_ > 1) //skip the first 2 layers
      	    {
      	      Tracking::Track aTrack=track;

      	      Tracking::TrackMeasurement aHit(*it);
	      aHit.hodo_=&alignedHodo;

	      //refit track without this hit and get residual at this layer. Track use still the "original hodo"
	      aTrack.fitAngle_=false;
	      aTrack.removeMeasurement(aTrack.hits_.begin()+iHit);
	      aTrack.fitTrack(); //now refit without the hit

      	      chi2+=aTrack.residual(aHit); //take residual at this layer
      	    }
	  ++iHit;
    	}
    }

  std::cout << "<RESIDUAL/TRACK>: " << chi2/tracks_.tracks_.size() << std::endl;
  return chi2/100.; // avoid fit instability
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
  int nFitParameters = (hodo_->layers_.size()-2)*4;
  for (int i=0;i<(hodo_->layers_.size()-2);++i) //keep the first 2 layers as reference
    {
      GlobalCoord_t layerPos=hodo_->layers_[i+2].position_;
      minimizer->SetLimitedVariable(i*4, Form("X_%d",i+2), layerPos(0), 1E-6, -30, 30);
      minimizer->SetLimitedVariable(i*4+1, Form("Y_%d",i+2), layerPos(1), 1E-6, -30, 30);
      if (alignZ_)
	minimizer->SetLimitedVariable(i*4+2, Form("Z_%d",i+2), layerPos(2), 1E-6, 0, 3000);
      else
	minimizer->SetFixedVariable(i*4+2, Form("Z_%d",i+2), layerPos(2));
      minimizer->SetLimitedVariable(i*4+3, Form("zRot_%d",i+2), hodo_->layers_[i+2].zRotation_, 1E-6, -0.5, 0.5);
    }


  //---fit
  ROOT::Math::Functor chi2(this, &AlignHodo::globalChi2, nFitParameters);
  minimizer->SetFunction(chi2);
  std::cout << "START ALIGNMENT MINIMIZATION WITH #" << tracks_.tracks_.size() << " TRACKS " << std::endl;
  clock_t begin = clock();
  minimizer->Minimize();
  clock_t end = clock();

  std::cout << "ELAPSED TIME: " << (double)(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;

  //---save aligned hodo
  for (int i=0;i<(hodo_->layers_.size()-2);++i) //keep the first 2 layers as reference
    {
      //--- set position and rotations
      for (int icoord=0;icoord<3;++icoord) 
	hodo_->layers_[i+2].position_(icoord)=minimizer->X()[i*4+icoord];
      hodo_->layers_[i+2].setZRotation(minimizer->X()[i*4+3]);
    }

  // hodo_->Print();

  // //---get covariance matrix
  // covarianceMatrixStatus_ = minimizer->CovMatrixStatus();
  // double covariances[nFitParameters * nFitParameters];
  // minimizer->GetCovMatrix(covariances);
  
  delete minimizer;        
}

bool AlignHodo::EndLoop(int iLoop, CfgManager& opts)
{
  minimize();
  return true;
}


bool AlignHodo::End(CfgManager& opts)
{
  // hodo_->SaveAs("align.root");
  return true;
}
