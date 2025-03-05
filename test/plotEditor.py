#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histfilename = "Xs_wdRcheck.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histoK")
histo.SetDirectory(0)
histfile.Close()


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.12)
canvas.SetRightMargin(0.08)
#canvas.SetTopMargin(0.2)
#canvas.SetLogy(True)


histo.SetAxisRange(3.5,6., "X")
#histo.SetAxisRange(1900, 8.e3, "Y")
histo.SetTitle("#it{J/#psi} + #it{K}^{#pm}; #it{M}_{inv} (GeV); # Events")
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

'''
funcfilename = "NPifuncs.root"
funcfile = r.TFile.Open(funcfilename)
Bfunc = r.gROOT.FindObject("Bfunc")
Xfunc = r.gROOT.FindObject("Xfunc")

canvas.cd()
Bfunc.SetLineColor(3)
Bfunc.Draw("same")
Xfunc.SetLineColor(2)
Xfunc.Draw("same")
'''

# Background line
bgdfilename = "JKdRcheckbgdfunc.root"
bgdfile = r.TFile.Open(bgdfilename)
bgd = r.gROOT.FindObject("bgd")
bgdparams = bgd.GetParameters()
expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2))+[3]"
bgdDrawFunc = r.TF1("bgdDrawFunc",expression,3.9,6.0,4)
bgdDrawFunc.SetParameters(bgdparams)

canvas.cd()
bgdDrawFunc.SetLineColor(2)
bgdDrawFunc.SetLineStyle(2)
bgdDrawFunc.Draw("same")

'''
l = r.TLatex()
l.SetTextFont(42)
l.SetTextSize(0.06)
l.DrawLatex(4.45,1100., "#it{X}_{#it{#pi}}^{#pm}")
l.DrawLatex(5.15,550.,"#it{B}^{#pm}")
l.SetTextSize(0.035)
l.DrawLatex(5.36,1180.,"|#it{M_{#mu#mu} - m_{J/#psi}}| < 0.1 GeV")
'''


canvas.Print("JKdRcheckbgd.pdf")
input('press enter to exit')