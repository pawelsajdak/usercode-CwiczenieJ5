#!/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_14_0_2/external/el9_amd64_gcc12/bin/python3
import ROOT as r
import sys

peakname = "18IV"
xmin = 0.01
xmax = 0.2
#par0 = 100.
#axmin = 3.5
#axmax = 6.

##########################################
histfilename = "BpmLTh.root"
histfile = r.TFile.Open(histfilename,"READ")
histo = histfile.Get("histo")
histo.SetDirectory(0)
histfile.Close()

#expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]*exp((-(x-[4])**2)/(2*[5]**2))+[6]"
fitFunc = r.TF1("fitFunc","expo",xmin,xmax)
fitFunc.SetParameters(6.5,-28.5)

results = histo.Fit(fitFunc,"ERSLB")
#funcFile = r.TFile.Open("NPifuncs.root","UPDATE")
#fitFunc.Write("Bfunc")
#funcFile.Close()

#with open('NPiresults.txt','a') as of:
    #print(peakname,"\t",fitFunc.GetParameter(1),"\n", results, file=of)

lifetime = (-1/fitFunc.GetParameter(1))/3.e10
print("Lifetime: ",lifetime," +/- ",3.e10*fitFunc.GetParError(1)*lifetime*lifetime)

canvas = r.TCanvas("canvas")
canvas.cd()
#canvas.SetLogy(True)

#histo.SetAxisRange(0.0,0.02)
#histo.SetAxisRange(3.5, 6., "X")
#histo.SetAxisRange(2.e3, 7.e3, "Y")
histo.SetTitle("Lifetime of B^{#pm};t [cm/c];Counts")
#histo.SetStats(0)
histo.Draw("h")


canvas.Print("temp.pdf")
input('press enter to exit')