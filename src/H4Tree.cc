#include "interface/H4Tree.h"

void H4Tree::Init()
{
    //---fill map< <board, group, channel>, pointer to first sample>
    //tree_->GetEntry(0);
    NextEntry();
    unsigned int currentDigiBoard=-1, currentDigiGroup=-1, currentDigiChannel=-1;
    for(unsigned int iSample=0; iSample<nDigiSamples; ++iSample)
    {
        if(digiBoard[iSample] != -1 && 
           (digiChannel[iSample] != currentDigiChannel ||
            digiGroup[iSample] != currentDigiGroup ||
            digiBoard[iSample] != currentDigiBoard))
        {
            currentDigiChannel = digiChannel[iSample];
            currentDigiGroup = digiGroup[iSample];
            currentDigiBoard = digiBoard[iSample];
            digiMap[make_tuple(currentDigiBoard, currentDigiGroup, currentDigiChannel)] = iSample;
            digiNSamplesMap[make_tuple(currentDigiBoard, currentDigiGroup, currentDigiChannel)] = 1;
        }
        else
            digiNSamplesMap[make_tuple(currentDigiBoard, currentDigiGroup, currentDigiChannel)]++;
    }

    NextEntry(GetTTreePtr()->GetEntriesFast()+1);
}

H4Tree::~H4Tree()
{
    delete[] evtTimeBoard;
    delete[] evtTime;
    delete[] adcBoard;
    delete[] adcChannel;
    delete[] adcData;
    delete[] tdcChannel;
    delete[] tdcData;
    delete[] pattern;
    delete[] patternBoard;
    delete[] patternChannel;
    delete[] triggerWords;
    delete[] triggerWordsBoard;
    delete[] digiBoard;
    delete[] digiGroup;
    delete[] digiChannel;
    delete[] digiSampleValue;
    delete[] digiSampleGain;
}
