#include "AlignTelescope.h"
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
bool AlignTelescope::Begin(CfgManager& opts, uint64* index)
{  
    //---inputs---
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstance");
  
    tLayout_=NULL;
    tracks_.tracks_.clear();
    return true;
}

//----------Begin-------------------------------------------------------------------------
bool AlignTelescope::BeginLoop(int iLoop, CfgManager& opts)
{  
    tracks_.tracks_.clear();
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool AlignTelescope::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---load tracks from trackReco
    string objectName="tracks";
    auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_tracks", "", false);

    Tracking::TrackContainer* tracks;

    if(shared_data.size() != 0)
        tracks = (Tracking::TrackContainer*)shared_data.at(0).obj;
    else
        cout << "[AlignTelescope::" << instanceName_ << "]: " << srcInstance_+"_tracks not found" << endl; 

    //---select good tracks
    for (auto& track : tracks->tracks_)
    {
        if (!tLayout_)
            tLayout_=track.tLayout_;

        if (track.hits_.size()>3) 
            tracks_.tracks_.push_back(track);
    }
  
    return true;
}

//---Compute global alignment Chi2
double AlignTelescope::globalChi2(const double* par)
{
    double chi2=0;

    Tracking::TelescopeLayout alignedTelescope(*tLayout_); //reset telescope
    std::vector<Tracking::TelescopeLayer> originalLayers;
    std::vector<Tracking::TelescopeLayer> alignedLayers;
    originalLayers.assign(alignedTelescope.layers_.begin(),alignedTelescope.layers_.end());
    alignedLayers.assign(alignedTelescope.layers_.begin(),alignedTelescope.layers_.end());

    if (par) //--- set the actual fit parameters
    {
        for (int i=0;i<(tLayout_->layers_.size()-2);++i) //keep the first 2 layers as reference
	{
            //--- set position and rotations
            for (int icoord=0;icoord<3;++icoord) 
                alignedLayers[i+2].position_(icoord)=par[i*4+icoord];
            alignedLayers[i+2].setZRotation(par[i*4+3]);

            alignedTelescope.layers_[i+2]=alignedLayers[i+2]; //put misaligned layer in layout 
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

                Tracking::TrackHit aHit(*it);
                aHit.tLayout_=&alignedTelescope;

                //refit track without this hit and get residual at this layer. Track use still the "original telescope"
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

//---Minuit minimizer for global alignment fit
void AlignTelescope::minimize()
{
    //---setup minimization
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);
    minimizer->SetTolerance(1e-6);
    minimizer->SetPrintLevel(1);
    //---setup variables and ranges
    int nFitParameters = (tLayout_->layers_.size()-2)*4;
    for (int i=0;i<(tLayout_->layers_.size()-2);++i) //keep the first 2 layers as reference
    {
        GlobalCoord_t layerPos=tLayout_->layers_[i+2].position_;
        minimizer->SetLimitedVariable(i*4, Form("X_%d",i+2), layerPos(0), 1E-6, -30, 30);
        minimizer->SetLimitedVariable(i*4+1, Form("Y_%d",i+2), layerPos(1), 1E-6, -30, 30);
        minimizer->SetLimitedVariable(i*4+2, Form("Z_%d",i+2), layerPos(2), 1E-6, 0, 3000);
        minimizer->SetLimitedVariable(i*4+3, Form("zRot_%d",i+2), tLayout_->layers_[i+2].zRotation_, 1E-6, -0.5, 0.5);
    }


    //---fit
    ROOT::Math::Functor chi2(this, &AlignTelescope::globalChi2, nFitParameters);
    minimizer->SetFunction(chi2);
    std::cout << "START ALIGNMENT MINIMIZATION WITH #" << tracks_.tracks_.size() << " TRACKS " << std::endl;
    clock_t begin = clock();
    minimizer->Minimize();
    clock_t end = clock();

    std::cout << "ELAPSED TIME: " << (double)(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;

    //---save aligned telescope
    for (int i=0;i<(tLayout_->layers_.size()-2);++i) //keep the first 2 layers as reference
    {
        //--- set position and rotations
        for (int icoord=0;icoord<3;++icoord) 
            tLayout_->layers_[i+2].position_(icoord)=minimizer->X()[i*4+icoord];
        tLayout_->layers_[i+2].setZRotation(minimizer->X()[i*4+3]);
    }

    // tLayout_->Print();

    // //---get covariance matrix
    // covarianceMatrixStatus_ = minimizer->CovMatrixStatus();
    // double covariances[nFitParameters * nFitParameters];
    // minimizer->GetCovMatrix(covariances);
  
    delete minimizer;        
}

bool AlignTelescope::EndLoop(int iLoop, CfgManager& opts)
{
    minimize();
    return true;
}


bool AlignTelescope::End(CfgManager& opts)
{
    return true;
}
