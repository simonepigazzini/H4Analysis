#ifndef __TRG_TREE__
#define __TRG_TREE__

#include <string>
#include <vector>

#include "TTree.h"

using namespace std;

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************

/**
   The TrgTree contains the information on the trigger configuration for each
   event. Only one trigger bit is active for each event, with the choice either made by 
   the shifter at the start of the run or by the DAQ in case laser events were injected
   during the interspil (in this case a single run can contain events with different
   trigger type whereas in the first case all the events in a run share the same trigger
   configuration).

   The available trigger types are configured through TriggerTypeFilter and stored in 
   the tree as branches such that one can select events recorded with a specific trigger
   simply with `trg==PHYS` or `trg==LASER` for instace.
 */
class TrgTree
{
public: 
    //---ctors---
    TrgTree(){};
    TrgTree(uint64* idx, TTree* tree=NULL);
    //---dtor---
    ~TrgTree(){};

    //---utils---
    void Init(map<int, string>& triggers);
    void Fill() {tree_->Fill();};
    void Reset();
    
    TTree* tree_; 

    uint64*       index;
    vector<int>   trgs_;
    /**
       Active trigger bit for this event.
     */
    unsigned int  trg;
};

#endif
