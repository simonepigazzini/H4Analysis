import ROOT
from parser_utils import *
from array import array

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'Save Geometry into a separate root file')
    parser.add_argument('-i', '--input' , dest='input', default='', help='input file')
    parser.add_argument('-o', '--output' , dest='output', default='', help='output file')
    parser.add_argument('-n', '--iterations' , dest='iterations', default='', help='input object name')
    parser.add_argument('-c', '--channels' , dest='channel', default='', help='input object name')
    parser.add_argument('-w', '--write' , dest='write', default='', help='output object name')
    args = parser.parse_args ()

    ROOT.gSystem.Load("lib/libH4Analysis.so")

    _file0=ROOT.TFile.Open(args.input)

    f=ROOT.TFile.Open(args.output,"RECREATE")
    for gr in [0]:
        for ch in [ 0 ]:
            for i in range(0,int(args.iterations)):
                calib=_file0.Get("Iter_%d/CalibDigi_digiCalib"%i)
                hDeltaT=ROOT.TH1F("deltaT_Gr%d_Ch%d_Iter%d"%(gr,ch,i),"deltaT_Gr%d_Ch%d_Iter%d"%(gr,ch,i),200,-0.5,0.5)
                hDeltaV=ROOT.TH1F("deltaV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),"deltaV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),800,-40,40)
                hSlopeV=ROOT.TH1F("slopeV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),"slopeV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),400,-0.01,0.01)
                hQuadraticV=ROOT.TH1F("quadraticV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),"quadraticV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),400,-0.01,0.01)
                hSlopeV_DeltaV=ROOT.TH2F("slopeVdeltaV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),"slopeVdeltaV_Gr%d_Ch%d_Iter%d"%(gr,ch,i),400,-10,10,400,-0.01,0.01)
                x,deltaT, deltaV, slopeV, quadraticV = array('d'), array('d'), array('d'), array('d'), array('d')
                for icell, cal in enumerate(calib.channelCalibrations_[(gr,ch)].calibrations_):
                    x.append(icell)
                    deltaT.append(cal.deltaT)
                    deltaV.append(cal.deltaV)
                    slopeV.append(cal.slopeV)
                    quadraticV.append(cal.quadraticV)
                    hDeltaT.Fill(cal.deltaT)
                    hDeltaV.Fill(cal.deltaV)
                    hSlopeV.Fill(cal.slopeV)
                    hQuadraticV.Fill(cal.quadraticV)
                    hSlopeV_DeltaV.Fill(cal.deltaV,cal.slopeV)

                grDeltaT = ROOT.TGraph( len(x), x, deltaT )
                grDeltaV = ROOT.TGraph( len(x), x, deltaV )
                grSlopeV = ROOT.TGraph( len(x), x, slopeV )
                grQuadraticV = ROOT.TGraph( len(x), x, quadraticV )
                f.cd()        
                grDeltaT.Write("deltaT_Graph_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                grDeltaV.Write("deltaV_Graph_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                grSlopeV.Write("slopeV_Graph_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                grQuadraticV.Write("quadraticV_Graph_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                hDeltaT.Write("deltaT_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                hDeltaV.Write("deltaV_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                hSlopeV.Write("slopeV_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                hQuadraticV.Write("quadraticV_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                hSlopeV_DeltaV.Write("slopeVdeltaV_Gr%d_Ch%d_Iter%d"%(gr,ch,i))
                if (i==int(args.iterations)-1 and args.write != ""):
                    calib.Write(args.write)
                    fExample=ROOT.TF1("fExample","[0]+[1]*x+[2]*x*x",0,4096)
                    fExample.SetParameter(0,deltaV[10])
                    fExample.SetParameter(1,slopeV[10])
                    fExample.SetParameter(2,quadraticV[10])
                    fExample.SetLineColor(ROOT.kRed)
                    fExample.Write()
f.ls()
f.Close()


