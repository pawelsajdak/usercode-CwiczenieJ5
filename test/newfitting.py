#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "psi(4415)"
xmin = 4.1
xmax = 4.9
par0 = 10.e3
#axmin = 3.5
#axmax = 6.

##########################################
histfilename = "Bspectrum.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmin+xmax)/2,0.05,50.e3,1.,1.)

results = histo.Fit(fitFunc,"ERS")

with open('Bresults.txt','a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
histo.SetAxisRange(3.5, 6., "X")
histo.SetAxisRange(38.e3, 52.e3, "Y")
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")

canvas.Print("fit2_"+peakname+".pdf")
input('press enter to exit')