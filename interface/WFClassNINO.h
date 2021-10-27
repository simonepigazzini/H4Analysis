#ifndef __WFCLASS_NINO_H__
#define __WFCLASS_NINO_H__

#include "WFClass.h"

using namespace std;

class WFClassNINO : public WFClass
{
public:
    //---ctors---
    WFClassNINO() {};
    WFClassNINO(int polarity, float tUnit);

    //---getters---
    void                  AddSample(float sample) override;
    float                 GetSignalIntegral(int thr, int min) override;
};

#endif
