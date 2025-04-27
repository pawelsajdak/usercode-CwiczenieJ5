#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histname = "histoPi"
nbins = 10

##########################################
histfilename = "3partMinv.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get(histname)
histo.SetDirectory(0)
histfile.Close()

print("Old histo bins: ",histo.GetNbinsX())
rbhisto = histo.Rebin(nbins,"rbhisto")

#'''
outFile = r.TFile.Open("rebin"+histfilename,"UPDATE")
rbhisto.Write(histname+str(nbins))
outFile.Close()
#'''

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

rbhisto.SetAxisRange(4.,7.)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
#histo.SetTitle("Lifetime of B^{#pm};t [cm/c];Counts")
#histo.SetStats(0)
rbhisto.Draw("h")


canvas.Print("temp.pdf")
input('press enter to exit')