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
#canvas.SetLogy(True)
canvas.SetBottomMargin(0.12)
canvas.SetLeftMargin(0.12)


histo.SetAxisRange(8.8,10.8, "X")
histo.SetAxisRange(7.e3,18.e3, "Y")
histo.SetTitle("Dimuon Events; #it{M_{inv}} (GeV); # Events")
histo.SetStats(0)
r.gStyle.SetTitleFontSize(0.07)
histo.SetLabelSize(0.04,"XY")
#histo.SetNdivisions(212, "X")
histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetTitleOffset(1.0)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleOffset(1.25)
histo.GetXaxis().CenterTitle(True)
histo.SetFillColor(19)
histo.Draw("h")

###########
expression = "[0]*exp((-(x-[1])**2)/(2*[2]**2)) + [3]+x*[4]+x*x*[5]"
fitFunc = r.TF1("fitFunc",expression,9.,9.75)
fitFunc.SetParameters(5312., 9.44819,0.074437, -369268., 81052.8,-4327.05)
fitFunc.SetLineColor(4)
fitFunc.Draw("same")

upsi2Sfunc = r.TF1("upsi2Sfunc",expression,9.75,10.2)
upsi2Sfunc.SetParameters(1569.16, 10.009, 0.0776291, 9138.18, 313.348, -26.871  )
upsi2Sfunc.SetLineColor(8)
upsi2Sfunc.Draw("same")

upsi3Sfunc = r.TF1("upsi3Sfunc",expression,10.2,10.5)
upsi3Sfunc.SetParameters(882.235, 10.3462,  0.0664682,     375213., -68762.4, 3227.27)
upsi3Sfunc.SetLineColor(2)
upsi3Sfunc.Draw("same")




l = r.TLatex()
l.SetTextFont(42)

l.SetTextSize(0.05)
l.DrawLatex(9.4,16500.,"#varUpsilon(1S)")
l.DrawLatex(9.9,12000.,"#varUpsilon(2S)")
l.DrawLatex(10.3,11000.,"#varUpsilon(3S)")

canvas.Print("Upsilonhisto.pdf")
input('press enter to exit')