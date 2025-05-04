#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

hname = "JpsiLT"
xmin = 2.e-12
xmax = 8.e-12

##########################################
histfilename = "h"+hname+".root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

expression = "[0]*exp(-x/[1])"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax,2)
fitFunc.SetParameters(40000.,1.e-12)

results = histo.Fit(fitFunc,"ERSLB")
#'''
funcFile = r.TFile.Open(hname+"_fitFunc.root","UPDATE")
fitFunc.Write()
funcFile.Close()

with open(hname+'_fitResults.txt','a') as of:
    print(results, file=of)
#'''


canvas = r.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

#histo.SetAxisRange(0.0,0.02)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
histo.SetTitle("Lifetime of J/#psi;t [s];Counts")
#histo.SetStats(0)
histo.Draw("h")
fitFunc.Draw("same")


canvas.Print(hname+"_Fit.pdf")
input('press enter to exit')
