#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

class Background:
    def __call__(self, arr,par):
        if (fitRange.IsInside(arr[0])):        
            return par[0] + par[1]*arr[0]
        else:
            r.TF1.RejectPoint()
            return 0.0
##########################################
histfilename = "KK_B0s.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("hKaonKaon")
histo.SetDirectory(0)
histfile.Close()

outname = "KK_B0s_bgd"
xmin = 0.99
xmax = 1.08

fitRange = r.Fit.DataRange()
fitRange.AddRange(xmin,1.01)
fitRange.AddRange(1.04,xmax)
b = Background()
fitFunc = r.TF1("fitFunc",b,xmin,xmax,2)
#fitFunc.SetParameters(-4750.,0.32,6180.)

#r.Math.MinimizerOptions.SetDefaultTolerance(1.e-6)
results = histo.Fit(fitFunc,"ERSL")

funcFile = r.TFile.Open(outname+"Func.root","RECREATE")
fitFunc.Write("bgd")
#funcFile.Close()

with open(outname+'Fit.txt','a') as of:
    print(results, file=of)
    
canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
histo.SetAxisRange(xmin-0.1, xmax+0.1, "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
#histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print(outname+".pdf")
input('press enter to exit')