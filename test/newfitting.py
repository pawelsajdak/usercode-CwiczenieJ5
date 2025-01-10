#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "Jpsi"
xmin = 2.5
xmax = 3.4
par0 = 10.e5
axmin = 0.
axmax = 12.

##########################################
histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmin+xmax)/2,0.1,1.,1.,1.)
fitFunc.FixParameter(5,0.)
results = histo.Fit(fitFunc,"ERS")
'''
with open('results2.txt') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)
'''
canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

histo.SetAxisRange(axmin, axmax)
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")

canvas.Print("fi_"+peakname+".pdf")
input('press enter to exit')