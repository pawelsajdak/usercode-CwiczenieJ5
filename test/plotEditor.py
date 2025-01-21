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
#canvas.SetLogy(True)


histo.SetAxisRange(3.5, 6., "X")
histo.SetAxisRange(38.e3, 52.e3, "Y")

histo.SetTitle("J/psi + kaon (zoom 3.5-6); Minv (GeV); #events")
histo.SetStats(0)

histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetTitleOffset(0.8)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleOffset(1.)
histo.Draw("h")
l = r.TLatex()
l.SetTextSize(0.05)

l.DrawLatex(4.33,505.e2,"#psi(4415)")
l.DrawLatex(5.25,46.e3,"B^{+-}")

#histo.Draw("h")
canvas.Print("Bzoom3_6.pdf")
input('press enter to exit')