#!/bin/python3

import os
import subprocess

import ROOT

from ROOT import std
from parser_utils import *

vstring = std.vector(std.string)

def IterativeProfiling(wf):
    """
    Process 2D histogram to provide a smooth average pulse shape
    """

    prof = ROOT.TH1F("prof", "", wf.GetNbinsX(), wf.GetXaxis().GetXmin(), wf.GetXaxis().GetXmax())
    for ibin in range(1, wf.GetNbinsX()+1):
        h = ROOT.TH1F("tmp", "", wf.GetNbinsY(), wf.GetYaxis().GetXmin(), wf.GetYaxis().GetXmax())
        for ybin in range(1, wf.GetNbinsY()+1):
            h.SetBinContent(ybin, wf.GetBinContent(ibin,ybin))
            h.SetBinError(ybin, wf.GetBinError(ibin,ybin))

        deltaMean = 9999
        deltaRMS = 9999
        oldMean = h.GetMean()
        oldRMS = h.GetRMS()
        while(deltaMean>1E-4 or deltaRMS>1E-5):
            h.GetXaxis().SetRangeUser(oldMean-3*oldRMS, oldMean+3*oldRMS)
            newMean = h.GetMean()
            newRMS = h.GetRMS()
            deltaMean = abs(newMean-oldMean)
            deltaRMS = abs(newRMS-oldRMS)
            oldMean = newMean
            oldRMS = newRMS

        #prof.SetBinContent(ibin, max(h.GetMean(), 0.))
        prof.SetBinContent(ibin, h.GetMean())
        prof.SetBinError(ibin, h.GetMeanError())

        h.Delete()

    return prof

if __name__ == '__main__':
    """
    Run template generation
    """

    parser = argparse.ArgumentParser (description = 'Run pulse shape template generation using DFT template plugin')
    parser.add_argument('-r', '--runs' , action=customAction, help='run to be processed, either list or file')
    parser.add_argument('-c', '--channels' , action=customAction, help='channel names')
    parser.add_argument('-o', '--output' , default='templates_file.root', help='output file path')
    parser.add_argument('-m', '--max-events' , default='2000', help='maximum number of events')
    parser.add_argument('-t', '--template-cfg' , dest='template_cfg', default='../cfg/H4DAQ_base.cfg', help='Template cfg file')
    parser.add_argument('--cut', default='1', help='Selection on input waveform')
    parser.add_argument('--bins', action=customAction, default=[1000, 0., 0.], help='Template histogram time binning')
    parser.add_argument('--dftt-instance', dest='dftt_instance', default='DFTTmpl', help='Instance name of the DFTTemplate plugin')
    parser.add_argument('--dryrun' , action="store_true", default=False, help='do not submit the jobs, just create them')
    parser.add_argument('--debug' , action="store_true", default=False, help='print debug messages')    
    
    args = parser.parse_args ()

    cfg = ROOT.CfgManager(args.template_cfg)

    for ch, run in zip(args.channels, args.runs):
        print('>>> Generating cfg for channel', ch, '...')
        base_dir = os.path.abspath(__file__)[:os.path.abspath(__file__).find('H4Analysis/')+len('H4Analysis/')]
        cfg.SetOpt(args.dftt_instance+'.channelsNames', vstring(1, ch))
        cfg.SetOpt(args.dftt_instance+'.outWFSuffix', vstring(1, '_T'))
        if args.debug:            
            template_cfg.Print(args.dftt_instance)
        new_cfg_file = base_dir+'tmp/'+args.template_cfg.strip('.cfg')+'_'+ch+'_TMPL.cfg'
        cfg.WriteToFile(new_cfg_file, True)

        print('>>> Running reconstruction on run', run, '...')
        cmd = base_dir+'bin/H4Reco '+new_cfg_file+' '+run
        ret = subprocess.getstatusoutput(cmd)
        if args.debug:
            print("Cmd:", cmd)
            print("Return:", ret)

        print('>>> Building the template for ', ch, '...')
        reco_file = ROOT.TFile.Open(base_dir+"/"+cfg.GetOpt('h4reco.outNameSuffix')+run+'.root')
        reco_tree = reco_file.Get("h4")

        tmpl_file = ROOT.TFile.Open(args.output, 'UPDATE')
        tmpl_file.Delete('tmpl_'+ch+';*')
        #h_tmpl = ROOT.TProfile('tmpl_'+ch, '', int(args.bins[0]), float(args.bins[1]), float(args.bins[2]), -0.1, 1.1)
        h_wf_2d = ROOT.TH2D('wf_2d_'+ch, '', int(args.bins[0]), float(args.bins[1]), float(args.bins[2]), 10000, -0.1, 1.1)
        
        var = 'wf_t.WF_val/digi_t.amp_max[%s]:wf_t.WF_time-digi_t.time[%s+CFD]' % (ch+'_T', ch+'_T')
        cut = 'wf_t.WF_ch=='+ch+'_T'+' && '+args.cut
        entries = reco_tree.Draw(var+'>>wf_2d_'+ch, cut)
        
        h_tmpl = IterativeProfiling(h_wf_2d)
        h_tmpl.Scale(1/h_tmpl.GetMaximum())        
        h_tmpl.Write('tmpl_'+ch)
        tmpl_file.Close()
        reco_file.Close()
        
