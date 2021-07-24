#include "HodoReco.h"

//**********Utils*************************************************************************
//----------Begin*************************************************************************
bool HodoReco::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    hodoFiberOrderA_.clear();
    hodoFiberOrderB_.clear();
  
    hodoFiberOrderA_.push_back(31);
    hodoFiberOrderA_.push_back(29);
    hodoFiberOrderA_.push_back(23);
    hodoFiberOrderA_.push_back(21);
    hodoFiberOrderA_.push_back(5);
    hodoFiberOrderA_.push_back(7);
    hodoFiberOrderA_.push_back(15);
    hodoFiberOrderA_.push_back(13);
    hodoFiberOrderA_.push_back(1);
    hodoFiberOrderA_.push_back(3);
    hodoFiberOrderA_.push_back(11);
    hodoFiberOrderA_.push_back(9);
    hodoFiberOrderA_.push_back(6);
    hodoFiberOrderA_.push_back(8);
    hodoFiberOrderA_.push_back(16);
    hodoFiberOrderA_.push_back(14);
    hodoFiberOrderA_.push_back(17);
    hodoFiberOrderA_.push_back(27);
    hodoFiberOrderA_.push_back(19);
    hodoFiberOrderA_.push_back(25);
    hodoFiberOrderA_.push_back(24);
    hodoFiberOrderA_.push_back(22);
    hodoFiberOrderA_.push_back(32);
    hodoFiberOrderA_.push_back(30);
    hodoFiberOrderA_.push_back(4);
    hodoFiberOrderA_.push_back(2);
    hodoFiberOrderA_.push_back(12);
    hodoFiberOrderA_.push_back(10);
    hodoFiberOrderA_.push_back(20);
    hodoFiberOrderA_.push_back(18);
    hodoFiberOrderA_.push_back(28);
    hodoFiberOrderA_.push_back(26);

    hodoFiberOrderB_.push_back(54);
    hodoFiberOrderB_.push_back(64);
    hodoFiberOrderB_.push_back(56);
    hodoFiberOrderB_.push_back(62);
    hodoFiberOrderB_.push_back(49);
    hodoFiberOrderB_.push_back(51);
    hodoFiberOrderB_.push_back(59);
    hodoFiberOrderB_.push_back(57);
    hodoFiberOrderB_.push_back(53);
    hodoFiberOrderB_.push_back(55);
    hodoFiberOrderB_.push_back(63);
    hodoFiberOrderB_.push_back(61);
    hodoFiberOrderB_.push_back(45);
    hodoFiberOrderB_.push_back(47);
    hodoFiberOrderB_.push_back(37);
    hodoFiberOrderB_.push_back(39);
    hodoFiberOrderB_.push_back(34);
    hodoFiberOrderB_.push_back(42);
    hodoFiberOrderB_.push_back(36);
    hodoFiberOrderB_.push_back(44);
    hodoFiberOrderB_.push_back(50);
    hodoFiberOrderB_.push_back(52);
    hodoFiberOrderB_.push_back(58);
    hodoFiberOrderB_.push_back(60);
    hodoFiberOrderB_.push_back(38);
    hodoFiberOrderB_.push_back(48);
    hodoFiberOrderB_.push_back(40);
    hodoFiberOrderB_.push_back(46);
    hodoFiberOrderB_.push_back(41);
    hodoFiberOrderB_.push_back(43);
    hodoFiberOrderB_.push_back(33);
    hodoFiberOrderB_.push_back(35);

    //---create a position tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    //---Register output trees (one per X or Y layer) if specified
    RegisterSharedData(new TTree("h1X", "hodo1X_tree"), "hodo1X_tree", storeTree);
    hodoTrees_[0] = PositionTree(index, (TTree*)data_.back().obj);
    hodoTrees_[0].Init();
    RegisterSharedData(new TTree("h1Y", "hodo1Y_tree"), "hodo1Y_tree", storeTree);
    hodoTrees_[1] = PositionTree(index, (TTree*)data_.back().obj);
    hodoTrees_[1].Init();
    RegisterSharedData(new TTree("h2X", "hodo2X_tree"), "hodo2X_tree", storeTree);
    hodoTrees_[2] = PositionTree(index, (TTree*)data_.back().obj);
    hodoTrees_[2].Init();
    RegisterSharedData(new TTree("h2Y", "hodo2Y_tree"), "hodo2Y_tree", storeTree);
    hodoTrees_[3] = PositionTree(index, (TTree*)data_.back().obj);
    hodoTrees_[3].Init();

    //---Register volatile hodo hit for track reconstruction
    for (int i=0; i<nPlanes_; ++i)
        RegisterSharedData(&hodoHits_[i], "hodo_layer_"+to_string(i), false);
    
    //---cluster size options
    if(!opts.OptExist(instanceName_+".minClusterSize") || !opts.OptExist(instanceName_+".maxClusterSize"))
    {
        cout << ">>> HodoReco ERROR: please specify minClusterSize and maxClusterSize" << endl;
        return false;
    }
    else
    {
        minClusterSize_ = opts.GetOpt<int>(instanceName_+".minClusterSize");
        maxClusterSize_ = opts.GetOpt<int>(instanceName_+".maxClusterSize");
    }

    return true;
}

bool HodoReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---clear output tree
    for(int i=0; i<nPlanes_; ++i)
    {
        hodoTrees_[i].Clear();
        hodoHits_[i].hits_.clear();
    }

    std::map<int,std::map<int,bool> > hodoFiberOn;
    
    for(int i=0; i<nPlanes_; ++i)
        for(int j=0; j<nFibers_; ++j)
            hodoFiberOn[i][j] = 0;
  
    std::map<int,std::vector<int> > fibersOn;
  
    for(unsigned int i=0; i<h4Tree.nPatterns; ++i)
    {
        if(h4Tree.patternBoard[i] == 134348801 ||
           h4Tree.patternBoard[i] == 134348802)
        {
            int pos = -1; // here is where the real hodoscope mapping is done
      
            if(h4Tree.patternBoard[i] == 134348801)
                pos = (h4Tree.patternChannel[i]<2) ? HODO_Y2 : HODO_X2;
            else if(h4Tree.patternBoard[i] == 134348802)
                pos = (h4Tree.patternChannel[i]<2) ? HODO_Y1 : HODO_X1;
      
            std::vector<int>* fiberorder = (bool)(h4Tree.patternChannel[i]&0b1) ? &hodoFiberOrderB_ : &hodoFiberOrderA_;
      
            for(unsigned int j=0; j<32; ++j)
            {
                bool thisfibon = (h4Tree.pattern[i]>>j)&0b1;
                hodoFiberOn[pos][fiberorder->at(j)-1] = thisfibon;
                if(thisfibon) fibersOn[pos].push_back(fiberorder->at(j)-1-32);

            }
        }
    }
    
    for(int i=0; i<nPlanes_; ++i)
    {
        std::sort(fibersOn[i].begin(), fibersOn[i].end());
        vector<vector<int> > clusters;
        for(auto& fiber : fibersOn[i])
        {
            if(clusters.size() > 0 && fiber-clusters.back().back() < 2)
                clusters.back().push_back(fiber);
            else
                clusters.emplace_back(vector<int>({fiber}));
        }
                
	
        Tracking::TelescopeLayout fakeHodo();
		
        for(auto& cluster : clusters)
        {
            if(cluster.size() >= minClusterSize_  && cluster.size() <= maxClusterSize_)
            {
                float value = 0;
                for(auto& fiber : cluster)
                    value += fiber;
                value /= cluster.size();

                if(i%2 == 0)
                {
                    Tracking::TrackHit trackMeasure(0.5*value, 0);
                    trackMeasure.setVarianceX(0.15*0.15);
                    trackMeasure.setVarianceY(9999.);
                    trackMeasure.calculateInverseVariance();
                    hodoHits_[i].hits_.push_back(trackMeasure);
                    hodoTrees_[i].n_clusters++;
                    hodoTrees_[i].clusters.emplace_back(cluster.size(), 0.5*value, -999, cluster.size());
                }
                else
                {
                    Tracking::TrackHit trackMeasure(0., -0.5*value);
                    trackMeasure.setVarianceX(9999.);
                    trackMeasure.setVarianceY(0.15*0.15);
                    trackMeasure.calculateInverseVariance();
                    hodoHits_[i].hits_.push_back(trackMeasure);
                    hodoTrees_[i].n_clusters++;
                    hodoTrees_[i].clusters.emplace_back(cluster.size(), -999, -0.5*value, cluster.size());
                }
            }
        }
    }
    
    //---fill output tree and clean up
    for(int i=0; i<nPlanes_; ++i)
        hodoTrees_[i].Fill();
    
    return true;
}
