import ROOT as r
import sys

peakname = "psi(2S)"
xmin = 3.5
xmax = 4.0
par0 = 30.e3



histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax)
fitFunc.SetParameters(par0,(xmin+xmax)/2,0.1,1.,1.,1.)
results = histo.Fit(fitFunc,"ERS")
'''
with open('results.txt', 'a') as of:
    print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)
'''
canvas = r.TCanvas("canvas")
canvas.cd()
histo.Draw("h")
canvas.Print("mycanvas.pdf")
input('press enter to exit')