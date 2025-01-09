import ROOT as r
import sys

histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,2.5,3.5)
fitFunc.SetParameters(700000.,3.,0.1,30000.,1.,1.)
histo.Fit(fitFunc,"ER")

canvas = r.TCanvas("canvas")
canvas.cd()
histo.Draw("h")
canvas.Print("mycanvas.pdf")
input('press enter to exit')