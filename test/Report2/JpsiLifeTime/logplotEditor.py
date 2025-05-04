#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histname = "JpsiLT"

histfilename = "h"+histname+".root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.12)
canvas.SetRightMargin(0.08)
#canvas.SetTopMargin(0.2)
canvas.SetLogy(True)


#histo.SetAxisRange(4.2,5.5, "X")
histo.SetAxisRange(8., 2.e5, "Y")
histo.SetTitle("Lifetime of J/#psi; Lifetime (s); Counts")
histo.SetStats(0)

r.gStyle.SetTitleFontSize(0.06)
histo.SetLabelSize(0.04,"XY")
#histo.SetNdivisions(40206, "X")
histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetTitleOffset(1.0)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleOffset(1.2)
histo.GetXaxis().CenterTitle(True)
histo.SetFillColor(19)
#histo.SetLineColor(28)
histo.Draw("h")

#'''
funcfilename = histname+"_fitFunc.root"
funcfile = r.TFile.Open(funcfilename,"READ")
fitFunc = r.gROOT.FindObject("fitFunc")

canvas.cd()
#fitFunc.SetLineColor(3)
fitFunc.Draw("same")
#'''

l = r.TLatex()
l.SetTextFont(42)
l.SetTextSize(0.035)
#l.SetTextAlign(31)
l.DrawLatex(56.e-13,8.e4,"#cbarM_{#mu#mu} - m_{J/#psi}#cbar < 0.02 GeV")
l.DrawLatex(565.e-14,3.e4,"#DeltaR(p_{J/#psi})_{ track,rec} < 0.15")


canvas.Print(histname+"_log"+"_finePlot.pdf")
input('press enter to exit')