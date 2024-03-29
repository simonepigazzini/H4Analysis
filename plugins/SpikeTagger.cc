#include "SpikeTagger.h"

//----------Utils-------------------------------------------------------------------------
bool SpikeTagger::Begin(map<string, PluginBase*>& plugins, CfgManager& opts, uint64* index)
{
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        Log("no source plugin specified", ERR);
        return false;
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");

    if (opts.OptExist(instanceName_+".weightsLd")) {
        weightsLd_ = opts.GetOpt<vector<float>>(instanceName_+".weightsLd");
    } else {
        // Use weights from TDR
        weightsLd_.assign({1.48322, -2.20018, 1.89766, -0.683441});
    }

    //---calculate channels for swiss cross and A_9 variable
    for (const auto& channel : channelsNames_) {
        channelsNames3By3_[channel].emplace_back(channel);
        const auto ch_letter_ascii = (int)channel[0];
        const auto ch_number_ascii = (int)channel[1];
        if (ch_letter_ascii >= int('B')) {
            auto ch_left(channel);
            ch_left[0] = char(int(ch_left[0]) - 1);
            if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_left) != channelsNames_.end()) {
                channelsNamesSwissCross_[channel].emplace_back(ch_left);
                channelsNames3By3_[channel].emplace_back(ch_left);
            }
            if (ch_number_ascii <= int('4')) {
                auto ch_left_top(ch_left);
                ch_left_top[1] = char(int(ch_left_top[1]) + 1);
                if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_left_top) != channelsNames_.end()) {
                    channelsNames3By3_[channel].emplace_back(ch_left_top);
                }
            }
            if (ch_number_ascii >= int('2')) {
                auto ch_left_bottom(ch_left);
                ch_left_bottom[1] = char(int(ch_left_bottom[1]) - 1);
                if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_left_bottom) != channelsNames_.end()) {
                    channelsNames3By3_[channel].emplace_back(ch_left_bottom);
                }
            }
        }
        if (ch_letter_ascii <= int('D')) {
            auto ch_right(channel);
            ch_right[0] = char(int(ch_right[0]) + 1);
            if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_right) != channelsNames_.end()) {
                channelsNamesSwissCross_[channel].emplace_back(ch_right);
                channelsNames3By3_[channel].emplace_back(ch_right);
            }
            if (ch_number_ascii <= int('4')) {
                auto ch_right_top(ch_right);
                ch_right_top[1] = char(int(ch_right_top[1]) + 1);
                if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_right_top) != channelsNames_.end()) {
                    channelsNames3By3_[channel].emplace_back(ch_right_top);
                }
            }
            if (ch_number_ascii >= int('2')) {
                auto ch_right_bottom(ch_right);
                ch_right_bottom[1] = char(int(ch_right_bottom[1]) - 1);
                if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_right_bottom) != channelsNames_.end()) {
                    channelsNames3By3_[channel].emplace_back(ch_right_bottom);
                }
            }
        }
        if (ch_number_ascii <= int('4')) {
            auto ch_top(channel);
            ch_top[1] = char(int(ch_top[1]) + 1);
            if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_top) != channelsNames_.end()) {
                channelsNamesSwissCross_[channel].emplace_back(ch_top);
                channelsNames3By3_[channel].emplace_back(ch_top);
            }
        }
        if (ch_number_ascii >= int('2')) {
            auto ch_bottom(channel);
            ch_bottom[1] = char(int(ch_bottom[1]) - 1);
            if (std::find(channelsNames_.begin(), channelsNames_.end(), ch_bottom) != channelsNames_.end()) {
                channelsNamesSwissCross_[channel].emplace_back(ch_bottom);
                channelsNames3By3_[channel].emplace_back(ch_bottom);
            }
        }
    }
    
    //---outputs---    
    string spikesTreeName = opts.OptExist(instanceName_+".spikeTreeName") ?
        opts.GetOpt<string>(instanceName_+".spikeTreeName") : "spikes";
    const bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree(spikesTreeName.c_str(), "spikes_tree"), "spikes_tree", storeTree);
    spikesTree_ = SpikesTree(index, (TTree*)data_.back().obj);
    spikesTree_.Init(channelsNames_);
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
    {
        int nSamples = 0;
        if(opts.OptExist(instanceName_+".storeNSampleAfterMax") && opts.OptExist(instanceName_+".storeNSampleBeforeMax"))
            nSamples = channelsNames_.size() *
                (opts.GetOpt<int>(instanceName_+".storeNSampleAfterMax") + opts.GetOpt<int>(instanceName_+".storeNSampleBeforeMax") + 1);
        else
        {
            for(const auto& channel : channelsNames_)
                nSamples += opts.GetOpt<int>(channel+".nSamples");
        }
        
        string wfTreeName = opts.OptExist(instanceName_+".wfTreeName") ?
            opts.GetOpt<string>(instanceName_+".wfTreeName") : "spikes_wf";
        RegisterSharedData(new TTree(wfTreeName.c_str(), "spikes_wf_tree"), "spikes_wf_tree", true);
        outWFTree_ = WFTree(nSamples, index, (TTree*)data_.back().obj);
        outWFTree_.Init();
    }

    return true;
}

bool SpikeTagger::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---setup output event 
    int outCh=0;
    bool fillWFtree=false;
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
        fillWFtree = *spikesTree_.index % opts.GetOpt<int>(instanceName_+".WFtreePrescale") == 0;

    //---load WFs from source instance shared data
    for(const auto& channel : channelsNames_)
    {
        const auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            Log("channels samples not found check DigiReco step", WARN); 
    }

    //---compute reco variables
    for(const auto& channel : channelsNames_)
    {
        //---skip dead channels
        if(WFs_.find(channel) == WFs_.end())
        {
            ++outCh;
            continue;
        }
        // Require at least 10 samples in the signal window just like in WFAnalyzer (it must not be less than in WFAnalyzer since otherwise the WFAnalyzer has dropped the channel and the variables are not available))
        if (WFs_[channel]->GetNSample() < opts.GetOpt<int>(channel+".signalWin", 0) + 10) {
            ++outCh;
            continue;
        }

        const auto analyzedWF = WFs_[channel]->GetSamples();
        const auto first_sample_time = WFs_[channel]->GetTimes()->at(0); // this might not be 0
        const auto t_unit = WFs_[channel]->GetTUnit();
        const int max_sample = static_cast<int>(std::round((WFs_[channel]->GetTimeCF(1).time - first_sample_time) / t_unit));

        //---Look for undershoot after maximum
        const auto undershoot_window = opts.GetOpt<int>(instanceName_+".undershootFinderWindow");
        if (max_sample + undershoot_window < analyzedWF->size()) 
        {
            auto undershoot_sample = std::min_element(analyzedWF->begin() + max_sample,
                                                      analyzedWF->begin() + max_sample + undershoot_window);
            spikesTree_.undershoot[outCh] = *undershoot_sample;
            spikesTree_.t_undershoot_minus_t_sample_max[outCh] = (std::distance(analyzedWF->begin(), undershoot_sample) - max_sample) * t_unit;
        }
        else 
        {
            spikesTree_.undershoot[outCh] = 1e5;
            spikesTree_.t_undershoot_minus_t_sample_max[outCh] = 1e5;
        }

        //---compute A/A_matrix locking the phase
        if(spikesTree_.max_hit == -1 || WFs_[channelsNames_[spikesTree_.max_hit]]->GetAmpMax() < WFs_[channel]->GetAmpMax())
            spikesTree_.max_hit = outCh;
        float matrix_A_sum = 0.;
        for(const auto& other_ch : channelsNames_) {
            if (max_sample < WFs_[other_ch]->GetSamples()->size()) { // in case the other channel has less samples than the channel with the overall maximum
                matrix_A_sum += WFs_[other_ch]->GetSamples()->at(max_sample);
            }
        }

        spikesTree_.amp_sum_matrix[outCh] = matrix_A_sum;

        //---Compute swiss cross variable 1 - A_4/A_1
        const float sample_max = analyzedWF->at(max_sample);
        float swiss_cross_A4_sum = 0.;
        const auto& channelNamesSwissCross = channelsNamesSwissCross_[channel];
        for(const auto& swiss_cross_ch : channelNamesSwissCross) {
            if (max_sample < WFs_[swiss_cross_ch]->GetSamples()->size()) { // in case the swiss cross channel has less samples than the channel with the overall maximum
                swiss_cross_A4_sum += WFs_[swiss_cross_ch]->GetSamples()->at(max_sample);
            }
        }

        float swiss_cross = -1e5;
        if (sample_max > 0) {
            swiss_cross = 1 - swiss_cross_A4_sum / sample_max;
        }

        spikesTree_.n_swiss_cross_neighbours[outCh] = channelNamesSwissCross.size();
        spikesTree_.swiss_cross[outCh] = swiss_cross;

        //---Compute 3 by 3 amplitude sum variable
        float amp_sum_3by3 = 0.;
        const auto& channel_names_3by3 = channelsNames3By3_[channel];
        for(const auto& ch_3by3 : channel_names_3by3) {
            if (max_sample < WFs_[ch_3by3]->GetSamples()->size()) { // in case the 3by3 channel has less samples than the channel with the overall maximum
                amp_sum_3by3 += WFs_[ch_3by3]->GetSamples()->at(max_sample);
            }
        }

        spikesTree_.n_channels_3by3[outCh] = channel_names_3by3.size();
        spikesTree_.amp_sum_3by3[outCh] = amp_sum_3by3;

        //---Compute number of samples over thresholds of 25, 50, and 75%
        int n_minus = 0;
        int n_plus = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            const float thr_frac = 0.75 - i * 0.25;
            const float thr = thr_frac * sample_max;
            while (max_sample - n_minus >= 0 and analyzedWF->at(max_sample - n_minus) > thr) {
                ++n_minus;
            }
            while (max_sample + n_plus + 1 < analyzedWF->size() and analyzedWF->at(max_sample + n_plus + 1) > thr) {
                ++n_plus;
            }
            switch (i) {
            case 0:
                spikesTree_.n_samples_above_75perc_max[outCh] = n_minus + n_plus;
                spikesTree_.tot_75perc_max[outCh] = (n_minus + n_plus) * t_unit;
                break;
            case 1:
                spikesTree_.n_samples_above_50perc_max[outCh] = n_minus + n_plus;
                spikesTree_.tot_50perc_max[outCh] = (n_minus + n_plus) * t_unit;
                break;
            case 2:
                spikesTree_.n_samples_above_25perc_max[outCh] = n_minus + n_plus;
                spikesTree_.tot_25perc_max[outCh] = (n_minus + n_plus) * t_unit;
            }
        }

        //---Compute ratios of sample amplitudes around the maximum to the maximum
        for (int i = -3; i <= 3; ++i) {
            if (i == 0)
                continue;

            if (max_sample + i > 0 and max_sample + i < analyzedWF->size() and sample_max > 0) {
                const float ratio = analyzedWF->at(max_sample + i) / sample_max;

                switch (i) {
                case -3:
                    spikesTree_.sample_max_minus3_over_sample_max[outCh] = ratio;
                    break;
                case -2:
                    spikesTree_.sample_max_minus2_over_sample_max[outCh] = ratio;
                    break;
                case -1:
                    spikesTree_.sample_max_minus1_over_sample_max[outCh] = ratio;
                    break;
                case 1:
                    spikesTree_.sample_max_plus1_over_sample_max[outCh] = ratio;
                    break;
                case 2:
                    spikesTree_.sample_max_plus2_over_sample_max[outCh] = ratio;
                    break;
                case 3:
                    spikesTree_.sample_max_plus3_over_sample_max[outCh] = ratio;
                }
            }
        }

        //---Calculate the linear discriminant LD value from the TDR (Eq. 9.2)
        const auto rminus1 = spikesTree_.sample_max_minus1_over_sample_max[outCh];
        auto ld = spikesTree_.sample_max_plus1_over_sample_max[outCh];
        float rminus1pow = 1.;
        for (const auto weightLd : weightsLd_) {
            ld -= weightLd * rminus1pow;
            rminus1pow *= rminus1;
        }
        spikesTree_.ld[outCh] = ld;

        //---Calculate delta t between maximum and 3 sigma of baseline rms
        const auto rms = WFs_[channel]->SubtractBaseline().rms;
        unsigned int s = 0;
        while (max_sample + s < analyzedWF->size() and analyzedWF->at(max_sample + s) > 3 * rms) {
            ++s;
        }
        spikesTree_.t_3sigma_noise_minus_t_sample_max[outCh] = s * t_unit;

        //---increase output tree channel counter
        ++outCh;        
    }


    //---WFs---
    if(fillWFtree)
    {
        outCh = 0;
        //---fix WF window around maximum sample of largest hit in the event.
        //   also assuming all channels have the same sampling frequency (i.e. don't mix VFEs with V1742 channels)
        const float tUnit = WFs_[channelsNames_[spikesTree_.max_hit]]->GetTUnit();
        const int max_sample = static_cast<int>(std::round(WFs_[channelsNames_[spikesTree_.max_hit]]->GetTimeCF(1).time / tUnit));
        for(const auto& channel : channelsNames_)
        {
            //---skip dead channels
            if(WFs_.find(channel) == WFs_.end())
            {
                ++outCh;
                continue;
            }
            const auto analyzedWF = WFs_[channel]->GetSamples();
            const auto times = WFs_[channel]->GetTimes();
            unsigned int firstSample = 0;
            unsigned int lastSample = analyzedWF->size();
            if(opts.OptExist(instanceName_+".storeNSampleAfterMax") && opts.OptExist(instanceName_+".storeNSampleBeforeMax"))
            {
                if (max_sample - opts.GetOpt<int>(instanceName_+".storeNSampleBeforeMax") >= 0) {
                    firstSample = max_sample - opts.GetOpt<int>(instanceName_+".storeNSampleBeforeMax");
                }
                lastSample = max_sample + opts.GetOpt<int>(instanceName_+".storeNSampleAfterMax");
            }
            for(unsigned int jSample=firstSample; jSample<=lastSample; ++jSample)
            {
                outWFTree_.WF_ch.push_back(outCh);
                if (jSample >= 0 && jSample < analyzedWF->size())
                    outWFTree_.WF_val.push_back(analyzedWF->at(jSample));
                else
                    outWFTree_.WF_val.push_back(0.);
                if (jSample >= 0 && jSample < times->size())
                    outWFTree_.WF_time.push_back(times->at(jSample));
                else
                    outWFTree_.WF_time.push_back(0.);
            }
            //---increase output tree channel counter
            ++outCh;        
        }
    }

    //---fill the output trees and clean
    //---reco var
    spikesTree_.Fill();
    spikesTree_.max_hit = -1;
    
    //---WFs
    if(fillWFtree)
        outWFTree_.Fill();

    return true;
}
