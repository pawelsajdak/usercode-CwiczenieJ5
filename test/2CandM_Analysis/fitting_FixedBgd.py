#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "KK"
xmin = 1.012
xmax = 1.07
##########################################
histfilename = "twoCandMass.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("hKaonKaon")
histo.SetDirectory(0)
histfile.Close()

# Background function with fitted parameters
bgdfilename = peakname+"bgdFunc.root"
bgdfile = r.TFile.Open(bgdfilename)
bgd = r.gROOT.FindObject("bgd")

# Fitting function
expression = "[0]+[1]*log(x) + [2]*exp((-(x-[3])**2)/(2*[4]**2))"
fitFunc = r.TF1("fitFunc",expression,xmin,xmax,5)
fitFunc.SetParameters(1.,1.,2000.,1.02,0.005)
fitFunc.FixParameter(0,bgd.GetParameter(0))
fitFunc.FixParameter(1,bgd.GetParameter(1))
#fitFunc.FixParameter(2,bgd.GetParameter(2))

#'''
results = histo.Fit(fitFunc,"ERSLB")
#funcFile = r.TFile.Open("NPifuncs.root","UPDATE")
#fitFunc.Write("Bfunc")
#funcFile.Close()

with open(peakname+'FitResults.txt','a') as of:
    print(results, file=of)
#'''
    
canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

histo.SetAxisRange(0.95,1.15)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(1500, 3.e3, "Y")
#histo.SetTitle("Lifetime of B^{#pm};t;Counts")
#histo.SetStats(0)
histo.Draw("h")
fitFunc.Draw("same")


canvas.Print(peakname+"Fit.pdf")
input('press enter to exit')
