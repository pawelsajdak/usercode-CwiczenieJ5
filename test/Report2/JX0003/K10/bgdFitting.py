#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

histname = "histoK10"

class Background:
    def __call__(self, arr,par):
        if (fitRange.IsInside(arr[0])):        
            return par[0] + par[1]*r.TMath.Exp(par[2]*arr[0])
        else:
            r.TF1.RejectPoint()
            return 0.0
##########################################
histfilename = "rebinJXnewdR.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get(histname)
histo.SetDirectory(0)
histfile.Close()

xmin = 4.4
xmax = 4.9
axmin = 4.0
axmax = 5.5

fitRange = r.Fit.DataRange()
fitRange.AddRange(xmin,4.56)
fitRange.AddRange(4.74,xmax)
b = Background()
fitFunc = r.TF1("fitFunc",b,xmin,xmax,3)
fitFunc.SetParameters(1300.,700000.,-1.4)

results = histo.Fit(fitFunc,"ERSL")

outname = histname+"hill"+"_bgd"
#'''
funcFile = r.TFile.Open(outname+"Func.root","RECREATE")
fitFunc.Write("bgd")
funcFile.Close()

with open(outname+'Fit.txt','a') as of:
    print(results, file=of)
#'''    

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

histo.SetAxisRange(axmin,axmax)
#histo.SetAxisRange(xmin-0.1, xmax+0.1, "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
#histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print(outname+".pdf")
input('press enter to exit')