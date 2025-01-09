import ROOT as r
import sys

peakname = "J/psi"
xmin = 2.5
xmax = 3.5
par0 = 600.e3



histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmax+xmin)/2,0.1,-300000.,200000.,-35000.)
results = histo.Fit(fitFunc,"ERS")

with open('results.txt', 'a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

canvas = r.TCanvas("canvas")
canvas.cd()
histo.Draw("h")
canvas.Print("mycanvas.pdf")
input('press enter to exit')