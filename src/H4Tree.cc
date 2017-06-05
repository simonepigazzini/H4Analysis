#include "interface/H4Tree.h"

void H4Tree::Init()
{
    //---fill map< <board, group, channel>, pointer to first sample>
    tree_->GetEntry(0);
    unsigned int currentDigiBoard=-1, currentDigiGroup=-1, currentDigiChannel=-1;
    for(unsigned int iSample=0; iSample<nDigiSamples; ++iSample)
    {
        if(digiChannel[iSample] != currentDigiChannel ||
           digiGroup[iSample] != currentDigiGroup ||
           digiBoard[iSample] != currentDigiBoard)
        {
            currentDigiChannel = digiChannel[iSample];
            currentDigiGroup = digiGroup[iSample];
            currentDigiBoard = digiBoard[iSample];            
            digiMap[make_tuple(currentDigiBoard, currentDigiGroup, currentDigiChannel)] = iSample;
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
