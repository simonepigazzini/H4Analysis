#include "TrackReco.h"
#include <regex>

#include "TRandom3.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"

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

  //---maxchi2 option
  maxChi2_ = opts.OptExist(instanceName_+".maxChi2") ?
    opts.GetOpt<float>(instanceName_+".maxChi2") : 100; 

  //---maxchi2 cleaning
  cleaningChi2Cut_ = opts.OptExist(instanceName_+".cleaningChi2Cut") ?
    opts.GetOpt<float>(instanceName_+".cleaningChi2Cut") : 100; 
  
  //---inputs---
  std::vector<string> layers = opts.GetOpt<vector<string> >(instanceName_+".layers");
  
  for (auto& layer: layers)
    {
      int measurementType=opts.GetOpt<int>(layer+".measurementType");
      std::vector<double> position=opts.GetOpt<vector<double> >(layer+".position");
      if (position.size() != 3)
	std::cout << "ERROR: Expecting a vector of size 3 for the layer position" << std::endl;
      GlobalCoord_t layerPos;
      layerPos.SetElements(position.begin(),position.end());
      
      //  RotationMatrix_t layer2Rot;
      ROOT::Math::Rotation3D::Scalar zRotationAngle = opts.GetOpt<ROOT::Math::Rotation3D::Scalar>(layer+".zRotationAngle");
      ROOT::Math::RotationZ r_z(zRotationAngle);      
      ROOT::Math::Rotation3D rotation(r_z);
      std::vector<double> rot_components(9);
      rotation.GetComponents(rot_components.begin());
      RotationMatrix_t layerRot(rot_components.begin(),rot_components.end());

      Tracking::TrackLayer aLayer(layerPos,layerRot);
      aLayer.measurementType_=measurementType;
      hodo_.addLayer(aLayer);      
    }

  hitProducers_ = opts.GetOpt<vector<string> >(instanceName_+".hitProducers");
  
  // layerMeas.reserve(hitProducers_.size());

  RegisterSharedData(&tracks_,"tracks",false);

  hodo_.Print();

  return true;
}

void TrackReco::buildTracks()
{
  while(1) 
    {
      Tracking::Track aTrack(&hodo_);  //create an empty track     
      aTrack.fitAngle_=false; //do not fit angle at this building step
      for (int i=0; i<hitProducers_.size(); ++i)
	{
	  std::string hitLayer=hitProducers_[i];

	  if (hits_[hitLayer]->hits_.size() == 0 )
	    continue;

	  //--- just add the first hit of a layer if no measurement 
	  if ( aTrack.hits_.size()==0                                               || //empty track 
	       (sqrt(aTrack.trackParCov_(0,0))>20. && hodo_.layers_[i].measureX() ) || //lousy track in X 
	       (sqrt(aTrack.trackParCov_(1,1))>20. && hodo_.layers_[i].measureY() )    //lousy track in Y 
	       )
	    {
	      Tracking::TrackMeasurement* hit=&(*hits_[hitLayer]->hits_.begin());
	      Tracking::TrackMeasurement aHit(hit->localPosition_,hit->localPositionError_,&hodo_,i); //create a new hit attached to the right geometry
	      aTrack.addMeasurement(aHit);
	      aTrack.fitTrack();
	      hits_[hitLayer]->hits_.erase( hits_[hitLayer]->hits_.begin() );
	      continue;
	    }

	  //--- otherwise look for the most compatible hit
	  std::vector<Tracking::TrackMeasurement>::iterator bestHit;
	  double minChi2=99999999.;

	  for (auto it=hits_[hitLayer]->hits_.begin(); it != hits_[hitLayer]->hits_.end(); ++it)
	    {	      
	      Tracking::TrackMeasurement aHit(it->localPosition_,it->localPositionError_,&hodo_,i); //create a new hit attached to the right geometry
	      double res=aTrack.residual(aHit);
	      if (res<minChi2)
		{
		  bestHit=it;
		  minChi2=res;
		}
	    }

	  if (minChi2<maxChi2_)
	    {
	      Tracking::TrackMeasurement aHit(bestHit->localPosition_,bestHit->localPositionError_,&hodo_,i); 
	      aTrack.addMeasurement(aHit);
	      aTrack.fitTrack();
	      hits_[hitLayer]->hits_.erase( bestHit );
	    }
	} //--end of loop over layers

      if (aTrack.hits_.size()==0)
	break; //no more Hits
      else
	tracks_.tracks_.push_back(aTrack);
    }
}

void TrackReco::cleanTracks()
{
  for (auto track=tracks_.tracks_.begin(); track !=tracks_.tracks_.end(); /* NOTHING */)
    {
      //want to have at least an hit on X and Y. Check the error
      if (track->covarianceMatrixStatus_ != 3 || 
	  sqrt(track->trackParCov_(0,0))>20.  || 
	  sqrt(track->trackParCov_(1,1))>20.  ||
	  track->chi2() > cleaningChi2Cut_)
  	{
	  tracks_.tracks_.erase( track );
  	}
      else
	++track;
    }
}

//----------ProcessEvent------------------------------------------------------------------
bool TrackReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
  //---load from source instances shared data
  for(auto& hitLayer : hitProducers_)
    {
      std::regex dot_re("\\.");
      std::sregex_token_iterator tkIter(hitLayer.begin(),hitLayer.end(),dot_re,-1);
      std::sregex_token_iterator tkIterEnd;
      std::vector<string> tokens;
      tokens.assign(tkIter,tkIterEnd);
      if (tokens.size() != 2)
	cout << "[TrackReco::" << instanceName_ << "]: Wrong input name " << hitLayer << endl;

      auto shared_data = plugins[tokens[0]]->GetSharedData(tokens[0]+"_"+tokens[1], "", false);

      if(shared_data.size() != 0)
	hits_[hitLayer] = (Tracking::LayerMeasurements*)shared_data.at(0).obj;
      else
	cout << "[TrackReco::" << instanceName_ << "]: " << tokens[0]+"_"+tokens[1] << " not found" << endl; 
    }

  tracks_.tracks_.clear();

  //---track building step  
  buildTracks();

  //---track cleaning
  cleanTracks();

  //---final fitting
  for (auto& track : tracks_.tracks_)
    {
      if (track.nFreeParameters()>4) //fit angle only when there is fitpix
	track.fitAngle_=true;
      else
	track.fitAngle_=false;
      track.fitTrack();
    }

  //---fill output tree
  trackTree_->Clear();
  trackTree_->n_tracks=tracks_.tracks_.size();
  for (auto& aTrack: tracks_.tracks_)
    {
      TrackPar par;
      par.value.assign(aTrack.trackPar_.Array(),aTrack.trackPar_.Array()+4);
      par.covariance.assign(aTrack.trackParCov_.Array(),aTrack.trackParCov_.Array()+10);
      trackTree_->fitResult.push_back(par);
      trackTree_->fitStatus.push_back(aTrack.covarianceMatrixStatus_);
      trackTree_->trackPattern.push_back(aTrack.trackPattern_);
      trackTree_->trackHits.push_back(aTrack.hits_.size());
      trackTree_->trackChi2.push_back(aTrack.chi2());
    }
  trackTree_->Fill();
  return true;
}
