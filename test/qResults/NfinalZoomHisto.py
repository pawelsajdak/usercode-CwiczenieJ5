#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histfilename = "fullhistogram.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.1)
#canvas.SetLogy(True)

histo.SetAxisRange(0.2, 1.4, "X")
histo.SetAxisRange(0., 130.e3, "Y")
histo.SetTitle("; #it{M_{inv}} (GeV); # Events")
histo.SetStats(0)
r.gStyle.SetTitleFontSize(0.07)
histo.SetLabelSize(0.04,"XY")
histo.SetNdivisions(40206, "X")
histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetTitleOffset(1.0)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleOffset(0.9)
histo.GetXaxis().CenterTitle(True)
histo.SetFillColor(19)
histo.Draw("h")
l = r.TLatex()
l.SetTextSize(0.055)

l.DrawLatex(0.53,50.e3,"#eta")
l.DrawLatex(0.75,115.e3,"#rho,#omega")
l.DrawLatex(1.,120.e3,"#phi")

canvas.Print("NNfinalZoomHisto.pdf")
input('press enter to exit')