#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys
import numpy as np

class Background:
    def __call__(self, arr,par):
        if (fitRange.IsInside(arr[0])):
            return par[0]*np.exp((-(arr[0]-par[1])**2)/(2*par[2]**2))+par[3]
        else:
            r.TF1.RejectPoint()
            return 0.0
##########################################
histfilename = "JxCorrN.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histoK")
histo.SetDirectory(0)
histfile.Close()

fitRange = r.Fit.DataRange()
fitRange.AddRange(3.8,4.2)
fitRange.AddRange(4.65,5.0)
fitRange.AddRange(5.5,5.9)
b = Background()
fitFunc = r.TF1("fitFunc",b,3.8,5.9,4)
fitFunc.SetParameters(1.e6,1.,1.,350.)

#r.Math.MinimizerOptions.SetDefaultTolerance(1.e-6)
results = histo.Fit(fitFunc,"ERSL")
funcFile = r.TFile.Open("NKbgd.root","RECREATE")
fitFunc.Write("bgd")
#funcFile.Close()

with open('NKbgd.txt','a') as of:
    print(results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
#histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print("NK_background.pdf")
input('press enter to exit')