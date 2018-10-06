#include "plugins/MakeCovarianceMatrix.h"

//----------Utils-------------------------------------------------------------------------
bool MakeCovarianceMatrix::Begin(CfgManager& opts, uint64* index)
{
    digiInstance_ = opts.GetOpt<string>(instanceName_+".digiInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");

    for(auto& ch : channelsNames_)
    {
        auto firstSample = opts.GetOpt<int>(ch+".templateFit.fitWin", 1);
        auto lastSample = opts.GetOpt<int>(ch+".templateFit.fitWin", 2);
        sums_[ch].resize(lastSample-firstSample);
        sum2s_[ch].resize(lastSample-firstSample);
        mapCovariances_[ch] = TH2F(ch.c_str(), "", lastSample+firstSample+1, -firstSample-0.5, lastSample+0.5,
                                   lastSample+firstSample+1, -firstSample-0.5, lastSample+0.5);
        mapCovariances_[ch].SetAxisRange(-1, 1, "Z");
        mapCovariances_[ch].SetContour(10000);
        RegisterSharedData(&mapCovariances_[ch], ch, true);
    }

    return true;
}

bool MakeCovarianceMatrix::ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    for(auto& channel : channelsNames_)
    {        
        WFClass* wf = (WFClass*)plugins[digiInstance_]->GetSharedData(digiInstance_+"_"+channel, "", false).at(0).obj;
        auto samples = wf->GetSamples();
        auto maxSample = wf->GetTimeCF(1.).first/wf->GetTUnit();
        auto firstSample = maxSample-opts.GetOpt<int>(channel+".templateFit.fitWin", 1);
        auto lastSample = maxSample+opts.GetOpt<int>(channel+".templateFit.fitWin", 2);
        if(wf->GetAmpMax() < 50.)
            continue;
        for(int iSample=firstSample; iSample<=lastSample; ++iSample)
        {
            values_[channel][iSample-firstSample].push_back(samples->at(iSample));
            sums_[channel][iSample-firstSample] += samples->at(iSample);
            sum2s_[channel][iSample-firstSample] += samples->at(iSample)*samples->at(iSample);
        }
        //---FIXME
        ++events_;   
    }

    return true;
}

bool MakeCovarianceMatrix::End(CfgManager& opts)
{
    for(auto& channel : channelsNames_)
    {        
        map<int, float> mean, rms;
        //---compute mean and rms of each sample
        for(unsigned int iSample=0; iSample<sums_[channel].size(); ++iSample)
        {
            mean[iSample] = sums_[channel][iSample]/events_;
            rms[iSample] = sqrt(sum2s_[channel][iSample]/events_ - mean[iSample]*mean[iSample]);
            // hMeans.Fill(mean[iSample]);
            // hRMSs.Fill(rms[iSample]);
        }
        //---compute covariances
        for(unsigned int iX=0; iX<sums_[channel].size(); ++iX)
        {
            for(unsigned int iY=0; iY<iX; ++iY)
            {
                float rho=0;
                for(int iEvt=0; iEvt<events_; ++iEvt)
                    rho += values_[channel][iX][iEvt]*values_[channel][iY][iEvt]
                        - values_[channel][iX][iEvt]*mean[iY] - values_[channel][iY][iEvt]*mean[iX];
                rho += events_*mean[iX]*mean[iY];
                rho = rho/((events_-1)*rms[iX]*rms[iY]);
                mapCovariances_[channel].SetBinContent(iX+1, iY+1, rho);
                mapCovariances_[channel].SetBinContent(iY+1, iX+1, rho);
            }
        }
    }
    
    return true;
}
