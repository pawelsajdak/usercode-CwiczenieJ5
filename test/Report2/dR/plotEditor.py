#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histfilename = "rebin_dRdistribution.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histodR20")
histo.SetDirectory(0)
zoomhisto = histfile.Get("histodR1")
zoomhisto.SetDirectory(0)
histfile.Close()


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.1)
canvas.SetRightMargin(0.04)
canvas.SetTopMargin(0.11)
#canvas.SetLogy(True)


histo.SetAxisRange(0.0,0.004, "X")
histo.SetAxisRange(0., 100., "Y")
histo.SetTitle("#DeltaR for M_{J/#psiK^{#pm}} #in [3.8, 6.0] GeV; min(#DeltaR); Counts")
histo.SetStats(0)

r.gStyle.SetTitleFontSize(0.06)
histo.SetLabelSize(0.04,"XY")
histo.SetNdivisions(505, "X")
histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetTitleOffset(1.0)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleOffset(0.9)
histo.GetXaxis().CenterTitle(True)
histo.SetFillColor(19)
#histo.SetLineColor(28)
histo.Draw("h")

##########  BOX #############
box = r.TBox(0.0,0.0,0.0008,10.)
box.SetFillStyle(0)
box.SetLineStyle(2)
box.SetLineColor(12)
box.SetLineWidth(3)
box.Draw("lsame")


##################  ZOOM    ###############
zoompad = r.TPad("zoompad","",0.2,0.25,0.9,0.88)
zoompad.SetLogy(1)
zoompad.Draw()
zoompad.cd()

zoompad.SetBottomMargin(0.07)
zoompad.SetLeftMargin(0.07)
zoompad.SetRightMargin(0.1)
zoompad.SetTopMargin(0.08)

zoomhisto.SetAxisRange(0.0,0.0008,"X")
zoomhisto.SetAxisRange(0.5,2.e4,"Y")
zoomhisto.SetTitle("Zoomed plot")
zoomhisto.SetStats(0)

r.gStyle.SetTitleFontSize(0.06)
zoomhisto.SetLabelSize(0.06,"XY")
#zoomhisto.SetNdivisions(40206, "X")
zoomhisto.GetXaxis().SetTitle("")
zoomhisto.GetYaxis().SetTitle("")
zoomhisto.SetFillColor(19)

zoomhisto.Draw()


canvas.Print("temp.pdf")
input('press enter to exit')