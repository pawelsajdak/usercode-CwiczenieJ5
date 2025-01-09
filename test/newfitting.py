import ROOT as r
import sys

peakname = "Y(3S)"
xmin = 10.2
xmax = 10.5
par0 = 10.e3
axmin = 8.
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
results = histo.Fit(fitFunc,"ERS")
#'''
with open('results.txt', 'a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)
#'''
canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

histo.SetAxisRange(axmin, axmax)
histo.SetTitle(peakname+"\t {:.3f}".format(fitFunc.GetParameter(1))+"; Minv; #events")
histo.SetStats(0)
histo.Draw("h")

canvas.Print("fit_"+peakname+".pdf")
input('press enter to exit')