#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histfilename = "Bspectrum.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

'''
histo.SetAxisRange(0.2, 1.4, "X")
histo.SetAxisRange(0., 130.e3, "Y")
histo.SetTitle("Dimuon Events; Minv (GeV); #events")
histo.SetStats(0)
histo.GetXaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleSize(0.05)
#histo.Draw("h")
l = r.TLatex()
l.SetTextSize(0.05)

l.DrawLatex(0.53,50.e3,"#eta")
l.DrawLatex(0.75,115.e3,"#rho,#omega")
l.DrawLatex(1.,120.e3,"#phi")
'''
histo.Draw("h")
canvas.Print("Bspectrum.pdf")
input('press enter to exit')