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
canvas.SetLogy(True)
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.1)


#histo.SetAxisRange(axmin, axmax)
histo.SetTitle("Dimuon Events; M_{inv} (GeV); #events")
histo.SetStats(0)
r.gStyle.SetTitleFontSize(0.07)
histo.SetLabelSize(0.04,"XY")
histo.SetNdivisions(212, "X")
histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetTitleOffset(1.0)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleOffset(0.9)
histo.GetXaxis().CenterTitle(True)
histo.SetFillColor(19)
histo.Draw("h")
l = r.TLatex()
l.SetTextSize(0.05)
l.SetTextFont(42)

l.DrawLatex(2.75,800.e3,"J/#psi")
l.DrawLatex(3.3,45.e3,"#psi(2S)")
l.DrawLatex(9.,20000.,"#varUpsilon(1,2,3S)")
l.DrawLatex(0.25,150.e3,"#eta #rho,#omega #phi")

canvas.Print("Nfinalhisto.pdf")
input('press enter to exit')