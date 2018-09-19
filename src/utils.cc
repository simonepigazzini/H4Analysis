#include "interface/utils.h"

//----------Smart Profile-----------------------------------------------------------------
void GetIterativeProfile(TH2* h2d, TH1* prof)
{
    //---loop over 2d-histo bins
    for (int ibin=1; ibin<=h2d->GetNbinsX(); ++ibin)
    {
        TH1F h_tmp("tmp", "", h2d->GetNbinsY(), h2d->GetYaxis()->GetXmin(), h2d->GetYaxis()->GetXmax());
        for (int ybin = 1;ybin<h2d->GetNbinsY()+1;++ybin)
        {
            h_tmp.SetBinContent(ybin, h2d->GetBinContent(ibin, ybin));
            h_tmp.SetBinError(ybin, h2d->GetBinError(ibin, ybin));
        }
        float deltaMean = 9999;
        float deltaRMS = 9999;
        float oldMean = h_tmp.GetMean();
        float oldRMS = h_tmp.GetRMS();
        while(deltaMean>1E-3 || deltaRMS>1E-4)
        {
            h_tmp.GetXaxis()->SetRangeUser(oldMean-3*oldRMS, oldMean+3*oldRMS);
            float newMean = h_tmp.GetMean();
            float newRMS = h_tmp.GetRMS();
            deltaMean = abs(newMean-oldMean);
            deltaRMS = abs(newRMS-oldRMS);
            oldMean = newMean;
            oldRMS = newRMS;
        }
        prof->SetBinContent(ibin, h_tmp.GetMean());
        prof->SetBinError(ibin, h_tmp.GetMeanError());

        h_tmp.Delete();
    }
    
    return;
}

