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

#histo.SetAxisRange(axmin, axmax)
histo.SetTitle("Dimuon Events; Minv (GeV); #events")
histo.SetStats(0)
histo.GetXaxis().SetTitleSize(0.05)
histo.GetYaxis().SetTitleSize(0.05)
histo.Draw("h")
l = r.TLatex()
l.SetTextSize(0.05)

l.DrawLatex(2.75,800.e3,"J/#psi")
l.DrawLatex(3.25,45.e3,"#psi(2S)")
l.DrawLatex(9.,20000.,"#varUpsilon(1,2,3S)")
l.DrawLatex(0.25,150.e3,"#eta #rho,#omega #phi")

canvas.Print("finalhisto.pdf")
input('press enter to exit')