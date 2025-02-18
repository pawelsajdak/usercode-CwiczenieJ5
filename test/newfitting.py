#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "XK"
xmin = 4.2
xmax = 4.8
par0 = 3000.
#axmin = 3.5
#axmax = 6.

##########################################
histfilename = "myVrtKPP.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histoK")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmin+xmax)/2,0.05,3.e3,1.,1.)

results = histo.Fit(fitFunc,"ERS")
funcFile = r.TFile.Open("Kfunctions.root","UPDATE")
fitFunc.Write("XK")
#funcFile.Close()

with open('Kresults.txt','a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#histo.SetAxisRange(axmin, axmax)
histo.SetAxisRange(3.5, 6., "X")
histo.SetAxisRange(2.e3, 8.e3, "Y")
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")


canvas.Print("KN_"+peakname+".pdf")
input('press enter to exit')