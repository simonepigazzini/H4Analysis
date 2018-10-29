#include "CalibDigitizer.h"
#include "TF1.h"
#include "TVectorD.h"

#include <time.h>

//----------Utils-------------------------------------------------------------------------
bool CalibDigitizer::Begin(CfgManager& opts, uint64* index)
{
    prevFit_=0.;
    thisFit_=0.;
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        cout << ">>> CalibDigitizer ERROR: no source plugin specified" << endl;
        return false;
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");

    nTotSamples_ = 0;
    for(auto& channel : channelsNames_)
        nTotSamples_ += opts.GetOpt<int>(channel+".nSamples");


    for(auto& channel : channelsNames_)
      {
	int nSamples=opts.GetOpt<int>(channel+".nSamples");
	pair<int,int> channelID(opts.GetOpt<int>(channel+".digiGroup"),opts.GetOpt<int>(channel+".digiChannel"));
	for (int iSample=0;iSample<nSamples;++iSample)
	  digiCalibration_.channelCalibrations_[channelID].calibrations_.push_back({0.,0.});
      }

    RegisterSharedData(&digiCalibration_, "digiCalib", true);

    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    string treeName = opts.OptExist(instanceName_+".treeName") ?
        opts.GetOpt<string>(instanceName_+".treeName") : "calib_tree";

    RegisterSharedData(new TTree(treeName.c_str(), treeName.c_str()), treeName.c_str(), storeTree);
  
    calibTree_ = new CalibDigiTree(index, (TTree*)data_.back().obj);
    calibTree_->Init();

    return true;
}


//----------Begin-------------------------------------------------------------------------
bool CalibDigitizer::BeginLoop(int iLoop, CfgManager& opts)
{  
  prevFit_=thisFit_;
  thisFit_=0;
  fitters_.clear();
  for (int i=0;i<nTotSamples_;++i)
    {
      fitters_.push_back(TLinearFitter(2)); //using a linear fitter with 2 parameters for each cell i % DeltaV_i + DeltaT_i*Slope %
      fitters_[i].SetFormula("hyp1"); 
    }
  
  return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool CalibDigitizer::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
  calibTree_->Clear();

    //---load WFs from source instance shared data
    for(auto& channel : channelsNames_)
    {
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            cout << "[CalibDigitizer::" << instanceName_ << "]: channels samples not found check DigiReco step" << endl; 
    }

    //---compute reco variables
    int offset=0;
    for(auto& channel : channelsNames_)
    {
        //---skip dead channels
        if(WFs_.find(channel) == WFs_.end())
            continue;

	pair<int,int> channelID(opts.GetOpt<int>(channel+".digiGroup"),opts.GetOpt<int>(channel+".digiChannel"));
	WFs_[channel]->SetCalibration(&digiCalibration_.channelCalibrations_[channelID]);
	WFs_[channel]->ApplyCalibration();
	
	auto analizedWF = WFs_[channel]->GetSamples();
	float tUnit = WFs_[channel]->GetTUnit();

	//Set function and fit initial values and ranges
	TF1* wFit=new TF1("wFit","[0]+[1]*TMath::Sin([2]*(x+[3]))",0,analizedWF->size()*tUnit);
	wFit->SetParameter(0,2400);
	wFit->SetParLimits(0,1800.,2600.);
	wFit->SetParameter(1,1350); 
	wFit->SetParLimits(1,900,1900);
	wFit->SetParameter(2,0.62); //80MHz 
	wFit->SetParLimits(2,0.5,0.7);
	wFit->SetParameter(3,0);
	wFit->SetParLimits(3,-50.,50.);

	double chi2=WFs_[channel]->AnalyticFit(wFit,0,analizedWF->size()-50); //fit using most of the samples avoiding the last samples of the digitizer

	if (chi2>50000)
	  continue;

	thisFit_+=chi2;

	//should add a check on the fit status
	//get residual and slope and add it to the linear minimizer
	const std::vector<double>* samples = WFs_[channel]->GetSamples();
	const std::vector<double>* times = WFs_[channel]->GetTimes();
	for (int iSample=0;iSample<samples->size();++iSample)
	  {
	    double slope=wFit->Derivative((*times)[iSample]*tUnit);
	    double residual=wFit->Eval((*times)[iSample]*tUnit)-(*samples)[iSample];
	    int cellIndex=(iSample + WFs_[channel]->GetStartIndexCell())%samples->size();
	    fitters_[offset+cellIndex].AddPoint(&slope,residual);
	  }	

	calibTree_->n_calibWf++;
	calibTree_->wfChi2.push_back(chi2);

	delete wFit;
	offset+=opts.GetOpt<int>(channel+".nSamples");
    }
  
    calibTree_->Fill();
    return true;
}

void CalibDigitizer::minimize()
{
  std::cout << "[CalibDigitizer]::Starting Mininimization for #" << nTotSamples_ << " cells" << std::endl;
  //  clock_t begin = clock();
  for (int i=0;i<nTotSamples_;++i)
    {
      //	fitters_[i].FixParameter(0,0.);
      fitters_[i].Eval();
    }
  //clock_t end = clock();
}

bool CalibDigitizer::EndLoop(int iLoop, CfgManager& opts)
{
  minimize();
  int offset=0;
  for(auto& channel : channelsNames_)
    {
      int nSamples=opts.GetOpt<int>(channel+".nSamples");
      auto tUnit = opts.GetOpt<float>(channel+".tUnit");
      pair<int,int> channelID(opts.GetOpt<int>(channel+".digiGroup"),opts.GetOpt<int>(channel+".digiChannel"));
      for (int iSample=0;iSample<nSamples;++iSample)
	{
	  TVectorD params;

	  fitters_[offset+iSample].GetParameters(params);
	  // TVectorD errors;
	  // fitters_[iSample].GetErrors(errors);

	  digiCalibration_.channelCalibrations_[channelID].calibrations_[iSample].deltaV+=params(0);
	  digiCalibration_.channelCalibrations_[channelID].calibrations_[iSample].deltaT+=(params(1)/tUnit);
	}
      offset+=nSamples;
    }

  std::cout << "[CalibDigitizer]::Relative Improvement " << fabs(thisFit_-prevFit_)/thisFit_ << std::endl;

  if (fabs(thisFit_-prevFit_)/thisFit_>1e-3)
    return true;
  else
    return false;
}


bool CalibDigitizer::End(CfgManager& opts)
{
    return true;
}

