import ROOT
from parser_utils import *
from array import array

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'Save Geometry into a separate root file')
    parser.add_argument('-i', '--input' , dest='input', default='', help='input file')
    parser.add_argument('-o', '--output' , dest='output', default='', help='output file')
    parser.add_argument('-n', '--iterations' , dest='iterations', default='', help='input object name')
    parser.add_argument('-c', '--channel' , dest='channel', default='', help='input object name')
#    parser.add_argument('-w', '--write' , dest='write', default='', help='output object name')
    args = parser.parse_args ()

    ROOT.gSystem.Load("lib/libH4Analysis.so")

    _file0=ROOT.TFile.Open(args.input)

    f=ROOT.TFile.Open(args.output,"RECREATE")    
    for i in range(0,int(args.iterations)):
        calib=_file0.Get("Iter_%d/CalibDigi_digiCalib"%i)
        hDeltaT=ROOT.TH1F("deltaT_Iter%d"%i,"deltaT_Iter%d"%i,200,-0.5,0.5)
        hDeltaV=ROOT.TH1F("deltaV_Iter%d"%i,"deltaV_Iter%d"%i,200,-2,2)
        x,deltaT, deltaV = array('d'), array('d'), array('d')
        for icell, cal in enumerate(calib.channelCalibrations_[(0,0)].calibrations_):
            x.append(icell)
            deltaT.append(cal.deltaT)
            deltaV.append(cal.deltaV)
            hDeltaT.Fill(cal.deltaT)
            hDeltaV.Fill(cal.deltaV)

        grDeltaT = ROOT.TGraph( len(x), x, deltaT )
        grDeltaV = ROOT.TGraph( len(x), x, deltaV )
        f.cd()        
        grDeltaT.Write("deltaT_Graph_Iter%d"%i)
        grDeltaV.Write("deltaV_Graph_Iter%d"%i)
        hDeltaT.Write("deltaT_Iter%d"%i)
        hDeltaV.Write("deltaV_Iter%d"%i)

f.ls()
f.Close()


