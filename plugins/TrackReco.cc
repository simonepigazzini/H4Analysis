#include "TrackReco.h"
#include <regex>

#include "TRandom3.h"

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
    opts.GetOpt<float>(instanceName_+".maxChi2") : 30; 
  
  //---inputs---
  std::vector<string> layers = opts.GetOpt<vector<string> >(instanceName_+".layers");
  
  for (auto& layer: layers)
    {
      int measurementType=opts.GetOpt<int>(layer+".measurementType");
      std::vector<double> position=opts.GetOpt<vector<double> >(layer+".position");
      //Need to take also rotation from cfg
      if (position.size() != 3)
	std::cout << "ERROR: Expecting a vector of size 3 for the layer position" << std::endl;
      GlobalCoord_t layerPos;
      layerPos.SetElements(position.begin(),position.end());
      Tracking::TrackLayer aLayer(layerPos);
      aLayer.measurementType_=measurementType;
      hodo_.addLayer(aLayer);

      //  RotationMatrix_t layer2Rot;
      //layer2Rot=ROOT::Math::SMatrixIdentity();
      // double angle=-0.15;
      // layer2Rot(0,0)=cos(angle);
      // layer2Rot(0,1)=sin(angle);
      // layer2Rot(1,0)=-sin(angle);
      // layer2Rot(1,1)=cos(angle);
    }

  hitProducers_ = opts.GetOpt<vector<string> >(instanceName_+".hitProducers");
  
  // layerMeas.reserve(hitProducers_.size());

  hodo_.Print();

  return true;
}

void TrackReco::buildTracks()
{
  while(1) 
    {
      std::cout << "+++++++++" << std::endl;
      Tracking::Track aTrack(&hodo_);  //create an empty track     
      aTrack.fitAngle_=false; //do not fit angle at this building step
      for (int i=0; i<hitProducers_.size(); ++i)
	{
	  std::string hitLayer=hitProducers_[i];
	  std::cout << "=======> " << i << " " <<  hits_[hitLayer]->hits_.size()  << std::endl;
	  for (auto it=hits_[hitLayer]->hits_.begin(); it != hits_[hitLayer]->hits_.end(); /* NOTHING */)
	    {
	      Tracking::TrackMeasurement* hit=&(*it);
	      Tracking::TrackMeasurement aHit(hit->localPosition_,hit->localPositionError_,&hodo_,i); //create a new hit attached to the right geometry
	      if ( aTrack.hits_.size()==0                                                        || //empty 
		   (sqrt(aTrack.trackParCov_(0,0))>20. && hodo_.layers_[i].measurementType_%2==1) || //lousy track in X 
		   (sqrt(aTrack.trackParCov_(1,1))>20. && hodo_.layers_[i].measurementType_%2==0) )  //lousy track in Y 
		{
		  aTrack.addMeasurement(aHit);
		  aTrack.fitTrack();
		  hits_[hitLayer]->hits_.erase( it );
		  break;
		}
	      else
		{
		  std::cout << "R " << aTrack.residual(aHit,false) << std::endl;
		  if (aTrack.residual(aHit)<maxChi2_)
		    {
		      aTrack.addMeasurement(aHit);
		      aTrack.fitTrack();
		      hits_[hitLayer]->hits_.erase( it );
		      break;
		    }
		  else
		    it++;
		}
	    }
	  if (aTrack.hits_.size()>0)
	      std::cout << sqrt(aTrack.trackParCov_(0,0)) << " " << sqrt(aTrack.trackParCov_(1,1)) << std::endl;
	}

      if (aTrack.hits_.size()==0)
	break; //no more Hits
      else
	{
	  tracks_.push_back(aTrack);
	  std::cout<< "OK " << aTrack.hits_.size() << std::endl;
	}
    }

}

void TrackReco::cleanTracks()
{
  for (auto track=tracks_.begin(); track !=tracks_.end(); /* NOTHING */)
    {
      //want to have at least an hit on X and Y. Check the error
      if (track->covarianceMatrixStatus_ != 3 || 
	  sqrt(track->trackParCov_(0,0))>20.  || 
	  sqrt(track->trackParCov_(1,1))>20.  ||
	  track->chi2() > 30.)
  	{
  	  std::cout << "ERASING TRACK with hits " << track->hits_.size() << std::endl;
	  tracks_.erase( track );
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
	cout << "[TrackReco::" << instanceName_ << "]: " << hitLayer << " not found" << endl; 
    }

  tracks_.clear();

  //---track building step  
  buildTracks();
  std::cout << "TRACKS AFTER BUILDING STEP " << tracks_.size() << std::endl;

  //---track cleaning
  cleanTracks();
  std::cout << "TRACKS AFTER CLEANING STEP " << tracks_.size() << std::endl;

  //---final fitting
  for (auto& track : tracks_)
    {
      if (track.hits_.size()>4)
	track.fitAngle_=true;
      else
	track.fitAngle_=false;
      track.fitTrack();
    }
  std::cout << "TRACKS AFTER FITTING STEP " << tracks_.size() << std::endl;

  //---fill output tree
  trackTree_->Clear();
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
