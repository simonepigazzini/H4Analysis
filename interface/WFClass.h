#ifndef __WFCLASS_H__
#define __WFCLASS_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;

#define MAX_INTERPOLATOR_POINTS 10000

struct WFBaseline
{
    float baseline;
    float rms;
    float slope;
    float k;
    float chi2;
};

struct WFFitResults
{
    double ampl;
    double time;
    double chi2;
    double slope;
};      

class  DigiChannelCalibration 
{
 public:

  struct calibrationParameters
  {
    double deltaV;
    double deltaT;
  };
  
  DigiChannelCalibration()
    {
      calibrations_.clear();
    }

  virtual ~DigiChannelCalibration()
    {
    }

  std::vector<calibrationParameters> calibrations_; //vector is indexed with startCellIndex

  ClassDef(DigiChannelCalibration, 1)
};

class DigitizerCalibration : public TNamed
{
 public:

 DigitizerCalibration() : TNamed()
    {
      channelCalibrations_.clear();
    }

  virtual ~DigitizerCalibration()
    {
    }
  
  map<pair<int,int>,DigiChannelCalibration> channelCalibrations_; 

  ClassDef(DigitizerCalibration, 1)
};


class WFClass : public TObject
{
public:
    //---ctors---
    WFClass() {};
    WFClass(int polarity, float tUnit);
    //---dtor---
    ~WFClass() {};

    //---getters---
    inline const vector<double>* GetSamples() {return &samples_;};
    inline int                   GetStartIndexCell() {return startIndexCell_;}
    inline int                   GetBWinMin() {return bWinMin_;}
    inline int                   GetBWinMax() {return bWinMax_;}
    inline int                   GetBIntWinMin() {return bIntWinMin_;}
    inline int                   GetBIntWinMax() {return bIntWinMax_;}
    inline float                 GetBaseline() {return baseline_;}
    inline float                 GetBaselineRMS() {return bRMS_;}
    inline int                   GetNSample() {return samples_.size();};
    inline int                   GetMaxSample() {return maxSample_;};
    inline float                 GetTUnit() {return tUnit_;};
    inline int                   GetSWinMin() {return sWinMin_;}
    inline int                   GetSWinMax() {return sWinMax_;}
    inline int                   GetSIntWinMin() {return sIntWinMin_;}
    inline int                   GetSIntWinMax() {return sIntWinMax_;}
    inline float                 GetFitAmpMax() {return fitAmpMax_;};
    inline float                 GetFitTimeMax() {return fitTimeMax_*tUnit_;};
    inline float                 GetLEThr() {return leThr_;};
    inline float                 GetLETime() {return leTime_;};
    inline float                 GetLESlope() {return leSlope_;};
    inline float                 GetTEThr() {return teThr_;};
    inline float                 GetTETime() {return teTime_;};
    inline float                 GetTESlope() {return teSlope_;};
    inline float                 GetCFFrac() {return cfFrac_;};
    inline float                 GetCFTime() {return cfTime_;};
    inline float                 GetCFSlope() {return cfSlope_;};
    inline TF1*                  GetAmpFunc() { return f_max_; };
    inline TF1*                  GetFitFunc() { return f_fit_; };
    float                        GetAmpMax(int min=-1, int max=-1);
    WFFitResults                 GetInterpolatedSample(int sample, int samplesLeft=-1, int samplesRight=-1);
    WFFitResults                 GetInterpolatedAmpMax(int min=-1, int max=-1, int nmFitSamples=7, int npFitSamples=7, string function="pol2", vector<float> params=vector<float>{});
    virtual WFFitResults         GetTime(string method, vector<float>& params); 
    virtual WFFitResults         GetTimeCF(float frac, int nFitSamples=5, int min=-1, int max=-1);
    virtual WFFitResults         GetTimeLE(float thr, int nmFitSamples=1, int npFitSamples=3, int min=-1, int max=-1);
    virtual WFFitResults         GetTimeTE(float thr, int nmFitSamples=1, int npFitSamples=3, int min=-1, int max=-1);
    float                        GetIntegral(int min=-1, int max=-1);
    float                        GetModIntegral(int min=-1, int max=-1);
    virtual float                GetSignalIntegral(int riseWin, int fallWin);
    //---setters---
    inline void                  SetTrigRef(float trigRef){trigRef_ = trigRef;};
    void                         SetSignalWindow(int min, int max);
    void                         SetSignalIntegralWindow(int min, int max);
    void                         SetBaselineWindow(int min, int max);
    void                         SetBaselineIntegralWindow(int min, int max);
    void                         SetTemplate(TH1* templateWF=NULL);
    void                         SetStartIndexCell(int cell){startIndexCell_=cell;};

    //---utils---
    void                         Reset();
    void                         AddSample(float sample) {samples_.push_back(polarity_*sample);};
    WFBaseline                   SubtractBaseline(int min=-1, int max=-1);
    WFFitResults                 TemplateFit(float offset=0., int lW=0, int hW=0);
    void                         AnalyticFit(TF1* f, int lW, int hW);
    void                         EmulatedWF(WFClass& wf, float rms, float amplitude, float time);
    void                         FFT(WFClass& wf, float tau, int cut);
    void                         Print();
    //---operators---
    WFClass&                     operator=(const WFClass& origin);
    WFClass                      operator-(const WFClass& sub);
    WFClass                      operator+(const WFClass& add);
    WFClass&                     operator-=(const WFClass& sub);
    WFClass&                     operator+=(const WFClass& add);
protected:
    //---utils---
    float                        BaselineRMS();
    float                        LinearInterpolation(float& A, float& B, const int& min, const int& max, const int& skipSample=-1);
    double                       TemplateChi2(const double* par=NULL);
    double                       AnalyticChi2(const double* par=NULL);
    
protected:
    vector<double> samples_;
    DigiChannelCalibration* calibration_;

    int            startIndexCell_;
    float          tUnit_;
    int            polarity_;
    float          trigRef_;
    int            sWinMin_;
    int            sWinMax_;
    int            sIntWinMin_;
    int            sIntWinMax_;
    int            bWinMin_;
    int            bWinMax_;
    int            bIntWinMin_;
    int            bIntWinMax_;
    int            maxSample_;
    float          fitAmpMax_;
    float          fitTimeMax_;
    float          fitChi2Max_;
    float          baseline_;
    float          bRMS_;
    int            cfSample_;
    float          cfFrac_;
    float          cfTime_;
    float          cfSlope_;
    int            leSample_;
    float          leThr_;
    float          leTime_;
    float          leSlope_;
    int            teSample_;
    float          teThr_;
    float          teTime_;
    float          teSlope_;
    float          chi2cf_;
    float          chi2le_;
    float          chi2te_;
    int            fWinMin_;
    int            fWinMax_;
    float          tempFitTime_;
    float          tempFitAmp_;
    TF1*           f_max_;
    TF1*           f_fit_;
    ROOT::Math::Interpolator* interpolator_;
};

#endif
