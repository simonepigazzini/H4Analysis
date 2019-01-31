#include "TrackReco.h"
#include <regex>

#include "TRandom3.h"
#include "Math/Rotation3D.h"
#include "Math/RotationZ.h"

//----------Begin-------------------------------------------------------------------------
bool TrackReco::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
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
    auto geoTag = opts.GetOpt<vector<string> >(instanceName_+".geometrySource");

    if(geoTag.size() > 1)
    {
        //---Load geometry from root file
        TFile* geoFile = TFile::Open(geoTag[0].c_str(), "READ");
        if(!geoFile)
        {
            cout << "[TrackReco::" << instanceName_ << "]: Cannot open file " << geoTag[0] << endl;
            return false;
        }            
        Tracking::TelescopeLayout* layout = (Tracking::TelescopeLayout*)geoFile->Get(geoTag[1].c_str());
        if(!layout)
        {
            cout << "[TrackReco::" << instanceName_ << "]: Cannot find object " << geoTag[1] << endl;
            return false;
        }

        tLayout_ = *((Tracking::TelescopeLayout*)layout->Clone("loadedLayout"));
        geoFile->Close();
    }
    else
    {
        //---load from cfg
        tLayout_ = Tracking::TelescopeLayout(opts, geoTag[0]);
    }

    hitProducers_ = opts.GetOpt<vector<string> >(instanceName_+".hitProducers");
  
    RegisterSharedData(&tracks_, "tracks", false);
    RegisterSharedData(&tLayout_, "telescope_layout", true);

    return true;
}

bool TrackReco::BeginLoop(int iLoop, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    tLayout_.Print();
    return true;
}

void TrackReco::buildTracks()
{
    while(1) 
    {
        Tracking::Track aTrack(&tLayout_);  //create an empty track     
        aTrack.fitAngle_=false; //do not fit angle at this building step
        for (int i=0; i<hitProducers_.size(); ++i)
        {
            std::string hitLayer=hitProducers_[i];
            if (hits_[hitLayer]->hits_.size() == 0 )
                continue;
            //--- just add the first hit of a layer if no measurement 
            if ( aTrack.hits_.size()==0                                               || //empty track 
                 (sqrt(aTrack.trackParCov_(0, 0))>20. && tLayout_.layers_[i].measureX() ) || //lousy track in X 
                 (sqrt(aTrack.trackParCov_(1, 1))>20. && tLayout_.layers_[i].measureY() )    //lousy track in Y 
                )
            {
                Tracking::TrackHit* hit=&(*hits_[hitLayer]->hits_.begin());
                Tracking::TrackHit aHit(hit->localPosition_, hit->localPositionError_, &tLayout_, i); //create a new hit attached to the right geometry
                aTrack.addMeasurement(aHit);
                aTrack.fitTrack();
                hits_[hitLayer]->hits_.erase( hits_[hitLayer]->hits_.begin() );
                continue;
            }
            //--- otherwise look for the most compatible hit
            std::vector<Tracking::TrackHit>::iterator bestHit;
            double minChi2=99999999.;
            for (auto it=hits_[hitLayer]->hits_.begin(); it != hits_[hitLayer]->hits_.end(); ++it)
            {	      
                Tracking::TrackHit aHit(it->localPosition_, it->localPositionError_, &tLayout_, i); //create a new hit attached to the right geometry
                double res=aTrack.residual(aHit);
                if (res<minChi2)
                {
                    bestHit=it;
                    minChi2=res;
                }
            }
            if (minChi2<maxChi2_)
            {
                Tracking::TrackHit aHit(bestHit->localPosition_, bestHit->localPositionError_, &tLayout_, i); 
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
    for(auto track = tracks_.tracks_.begin(); track != tracks_.tracks_.end(); /* NOTHING */)
    {
        //---want to have at least an hit on X and Y. Check the error
        if(track->covarianceMatrixStatus_ != 3 || 
           sqrt(track->trackParCov_(0, 0))>20.  || 
           sqrt(track->trackParCov_(1, 1))>20.  ||
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
        std::sregex_token_iterator tkIter(hitLayer.begin(), hitLayer.end(), dot_re, -1);
        std::sregex_token_iterator tkIterEnd;
        std::vector<string> tokens;
        tokens.assign(tkIter, tkIterEnd);
        if(tokens.size() != 2)
        {
            cout << "[TrackReco::" << instanceName_ << "]: Wrong input name " << hitLayer << endl;
            return false;
        }

        auto shared_data = plugins[tokens[0]]->GetSharedData(tokens[0]+"_"+tokens[1], "", false);

        if(shared_data.size() != 0)
            hits_[hitLayer] =(Tracking::LayerHits*)shared_data.at(0).obj;
        else
            cout << "[TrackReco::" << instanceName_ << "]: " << tokens[0]+"_"+tokens[1] << " not found" << endl; 
    }

    tracks_.tracks_.clear();

    //---track building step  
    buildTracks();

    //---track cleaning
    cleanTracks();
    //---final fitting
    for(auto& track : tracks_.tracks_)
    {
        if(track.nFreeParameters()>4) //fit angle only when there is fitpix
            track.fitAngle_=true;
        else
            track.fitAngle_=false;
        track.fitTrack();
    }

    //---fill output tree
    trackTree_->Clear();
    trackTree_->n_tracks=tracks_.tracks_.size();
    for(auto& aTrack: tracks_.tracks_)
    {
        TrackPar par;
        par.value.assign(aTrack.trackPar_.Array(), aTrack.trackPar_.Array()+4);
        par.covariance.assign(aTrack.trackParCov_.Array(), aTrack.trackParCov_.Array()+10);
        trackTree_->fitResult.push_back(par);
        trackTree_->fitStatus.push_back(aTrack.covarianceMatrixStatus_);
        trackTree_->trackPattern.push_back(aTrack.trackPattern_);
        trackTree_->trackHits.push_back(aTrack.hits_.size());
        trackTree_->trackChi2.push_back(aTrack.chi2());
    }

    trackTree_->Fill();
    return true;
}
