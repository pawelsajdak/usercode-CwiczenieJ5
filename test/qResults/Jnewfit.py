#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "Jpsi"
xmin = 2.2
xmax = 3.5
par0 = 1.e6-3.e4
#axmin = 3.5
#axmax = 6.

##########################################
histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]+x*x*x*[6]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
#expression = "[0] + x*[1]+x*x*[2]"
#fitFunc = r.TF1("fitFunc",expression,2.0,2.6)

fitFunc.SetParameters(par0,3.1,0.05,50.e3,1.,1.,1.)
fitFunc.FixParameter(3,42707.3)
fitFunc.FixParameter(4, 2236.93)
fitFunc.FixParameter(5, -1932.33)
fitFunc.FixParameter(6,  -402.025)

results = histo.Fit(fitFunc,"ERS")
funcFile = r.TFile.Open("Jfunction.root","RECREATE")
fitFunc.Write()
#funcFile.Close()

with open('Jresults.txt','a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
#histo.SetAxisRange(1.8, 2.65, "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print("NewFit"+peakname+".pdf")
input('press enter to exit')