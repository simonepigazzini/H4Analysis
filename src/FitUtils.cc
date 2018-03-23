#include "interface/FitUtils.h"



/*** double crystalBall ***/
double crystalBallLowHigh(double* x, double* par)
{
    //[0] = N
    //[1] = mean
    //[2] = sigma
    //[3] = alpha
    //[4] = n
    //[5] = alpha2
    //[6] = n2
  
    double xx = x[0];
    double mean = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double n = par[4];
    double alpha2 = par[5];
    double n2 = par[6];
  
    if( (xx-mean)/sigma > fabs(alpha) )
    {
        double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
        double B = n/fabs(alpha) - fabs(alpha);

        return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n);
    }

    else if( (xx-mean)/sigma < -1.*fabs(alpha2) )
    {
        double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
        double B = n2/fabs(alpha2) - fabs(alpha2);

        return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n2);
    }

    else
    {
        return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
    }
}



/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity)
{
    float integralMax = fraction * histo->Integral();
  
    int N = histo -> GetNbinsX();
    std::vector<float> binCenters(N);
    std::vector<float> binContents(N);
    std::vector<float> binIntegrals(N);
    for(int bin1 = 0; bin1 < N; ++bin1)
    {
        binCenters[bin1] = histo->GetBinCenter(bin1+1);
        binContents[bin1] = histo->GetBinContent(bin1+1);
    
        for(int bin2 = 0; bin2 <= bin1; ++bin2)
            binIntegrals[bin1] += binContents[bin2];
    }
  
    float min = 0.;
    float max = 0.;
    float delta = 999999.;
    for(int bin1 = 0; bin1 < N; ++bin1)
    {
        for(int bin2 = bin1+1; bin2 < N; ++bin2)
        {
            if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
            float tmpMin = histo -> GetBinCenter(bin1);
            float tmpMax = histo -> GetBinCenter(bin2);
      
            if( (tmpMax-tmpMin) < delta )
            {
                delta = (tmpMax - tmpMin);
                min = tmpMin;
                max = tmpMax;
            }
      
            break;
        }
    }
  
    TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
    for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
    {
        if( smallHisto->GetBinCenter(bin) < min )
            smallHisto -> SetBinContent(bin,0);
    
        if( smallHisto->GetBinCenter(bin) > max )
            smallHisto -> SetBinContent(bin,0);
    }
    smallHisto -> SetFillColor(kYellow);
  
    float mean = smallHisto -> GetMean();
    float meanErr = smallHisto -> GetMeanError();  
  
    ret[0] = mean;
    ret[1] = meanErr;
    ret[2] = min;
    ret[3] = max;
}
