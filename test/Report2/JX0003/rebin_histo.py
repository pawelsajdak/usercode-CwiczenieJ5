#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

histname = "histoPi"
nbins = 20

##########################################
histfilename = "JXnewdR.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get(histname)
histo.SetDirectory(0)
histfile.Close()

print("Old histo bins: ",histo.GetNbinsX())
rbhisto = histo.Rebin(nbins,"rbhisto")

rbhisto.SetAxisRange(4.,6.5)
#'''
outFile = r.TFile.Open("rebin"+histfilename,"UPDATE")
rbhisto.Write(histname+str(nbins))
outFile.Close()
#'''

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#rbhisto.SetAxisRange(4.,7.)

rbhisto.Draw("h")


canvas.Print("temp.pdf")
input('press enter to exit')