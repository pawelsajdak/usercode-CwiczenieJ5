#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histfilename = "myVrtKPP.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histoK")
histo.SetDirectory(0)
histfile.Close()


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.12)
#canvas.SetTopMargin(0.2)
#canvas.SetLogy(True)


histo.SetAxisRange(3.5,6., "X")
histo.SetAxisRange(1900, 8.e3, "Y")
histo.SetTitle("J/#psi + K^{+-}; #it{M_{inv}} (GeV); # Events")
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
histo.Draw("h")


funcfilename = "Kfunctions.root"
funcfile = r.TFile.Open(funcfilename)
Bfunc = r.gROOT.FindObject("fitFunc")
Xfunc = r.gROOT.FindObject("XK")

canvas.cd()
Bfunc.SetLineColor(3)
Bfunc.Draw("same")
Xfunc.SetLineColor(2)
Xfunc.Draw("same")
#canvas.Draw()

l = r.TLatex()
l.SetTextFont(42)
l.SetTextSize(0.06)
l.DrawLatex(4.42,6200., "X_{K}^{+-}")
l.DrawLatex(5.25,3100.,"B^{+-}")
l.SetTextSize(0.035)
l.DrawLatex(5.3,7400.,"#splitline{|M_{#mu#mu} - m_{J/#psi}| < 0.1 GeV}{6.5 #times 10^{6} entries}")



canvas.Print("KN.pdf")
input('press enter to exit')