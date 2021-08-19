#ifndef __WFCLASS_LITEDTU_H__
#define __WFCLASS_LITEDTU_H__

#include "WFClass.h"
#include "TGraph.h"

using namespace std;

class WFClassLiTEDTU : public WFClass
{
public:
    //---ctors---
    WFClassLiTEDTU() {};
    WFClassLiTEDTU(int polarity, float tUnit, DigiChannelCalibration* calibration=NULL);
    //---getters---
    void                           AddSample(float sample, float gain) override; 
    WFFitResults                   TemplateFit(float ampl_threshold=0, float offset=0., int lW=0, int hW=0) override;
};

#endif
