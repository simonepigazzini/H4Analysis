#include "interface/H4Tree.h"

void H4Tree::Init()
{
    //---fill map< <group, channel>, pointer to first sample>
    tree_->GetEntry(0);
    unsigned int currentDigiGroup=-1, currentDigiChannel=-1;
    for(unsigned int iSample=0; iSample<nDigiSamples; ++iSample)
    {
        if(digiChannel[iSample] != currentDigiChannel)
        {
            currentDigiChannel = digiChannel[iSample];
            if(digiGroup[iSample] != currentDigiGroup)
                currentDigiGroup = digiGroup[iSample];
            digiMap[make_pair(currentDigiGroup, currentDigiChannel)]=iSample;
        }
    }
}

H4Tree::~H4Tree()
{
    delete[] tdcChannel;
    delete[] tdcData;
    delete[] pattern;
    delete[] patternBoard;
    delete[] patternChannel;
    delete[] digiGroup;
    delete[] digiChannel;
    delete[] digiSampleValue;
}
