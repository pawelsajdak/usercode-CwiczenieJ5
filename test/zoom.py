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
histo.SetAxisRange(0.,1.5)

canvas = r.TCanvas("canvas")
canvas.cd()
histo.Draw("h")
canvas.Print("zoom.pdf")
input('press enter to exit')