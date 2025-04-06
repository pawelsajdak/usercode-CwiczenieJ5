#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

class Background:
    def __call__(self, arr,par):
        if (fitRange.IsInside(arr[0])):        
            return par[0] + par[2]*(arr[0]-par[1])**2
        else:
            r.TF1.RejectPoint()
            return 0.0
##########################################
histfilename = "twoCandMass.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("hKaonKaon")
histo.SetDirectory(0)
histfile.Close()

outname = "KKbgd"
xmin = 1.003
xmax = 1.1

fitRange = r.Fit.DataRange()
fitRange.AddRange(xmin,1.008)
fitRange.AddRange(1.035,xmax)
b = Background()
fitFunc = r.TF1("fitFunc",b,xmin,xmax,3)
fitFunc.SetParameters(5.e3,1.09,-500.e3)

#r.Math.MinimizerOptions.SetDefaultTolerance(1.e-6)
results = histo.Fit(fitFunc,"ERSL")

funcFile = r.TFile.Open(outname+"Func.root","RECREATE")
fitFunc.Write("bgd")
#funcFile.Close()

with open(outname+'Fit.txt','a') as of:
    print(results, file=of)
    
canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
histo.SetAxisRange(xmin, xmax, "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
#histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print(outname+".pdf")
input('press enter to exit')