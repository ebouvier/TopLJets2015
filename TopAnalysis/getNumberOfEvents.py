#! /usr/bin/env python

import math, os, string, sys

from ROOT import *

texDir = "sel"

plotDir = os.path.join(texDir, "plots")
if not os.path.isdir(texDir):
    parser.error("you must run the plotter before")

rootFile = TFile.Open(os.path.join(plotDir, "plotter.root"))    
procMC = [["\\ttbar", "t#bar{t}"], ["Single t", "t-ch"], ["\\PZ$+$jets", "DY"], ["\\PW$+$jets", "W"], ["t\\PW", "tW"], ["\\ttbar$+$V", "t#bar{t}+V"], ["$\\PW\\PW/\\PW\\PZ/\\PZ\\PZ$", "Multiboson"]]
hYields = [["in the e$+$jets channel", "yields_e_check"], ["in the $\\mu+$jets channel", "yields_m_check"], ["in the ee$+$jets channel", "yields_ee_check"], ["in the e$\\mu+$jets channel", "yields_em_check"], ["in the $\\mu\\mu+$jets", "yields_mm_check"], ["inclusively", "yields_all_check"]]
hDecays = [["\\ttbar\\", "decay_all_topsel"], ["$\\text{J}/\\psi$ meson", "decay_all_jpsicand"], ["$\\text{D}^{0}$ meson", "decay_all_d0cand"]]

texFile = "NumberOfEvents.tex"
tex = open(os.path.join(texDir,texFile), 'w')

tex.write("\\documentclass{article}\n")
tex.write("\\usepackage[a3paper, margin=10mm,landscape]{geometry}\n")
tex.write("\\usepackage{array}\n")
tex.write("\\usepackage{multirow}\n")
tex.write("\\usepackage{amsmath}\n")
tex.write("\\pagestyle{empty}\n")
tex.write("\\newcommand{\\ttbar}{\\ensuremath{\\text{t}\\bar{\\text{t}}}}\n")
tex.write("\\newcommand{\\PW}{\\ensuremath{\\text{W}}}\n")
tex.write("\\newcommand{\\PZ}{\\ensuremath{\\text{Z}}}\n")
tex.write("\\begin{document}\n")
tex.write("\n%%%\n")


tex.write("\n\n\\section{Step by step}")

tex.write("\n\\subsection{Event yields}\n")    

for hYield in hYields:
    tYield = hYield[0]
    nYield = hYield[1]
    tex.write("\\begin{table}[h!]\n")
    tex.write("\\begin{center}\n")
    tex.write("\\caption{\\label{%s}\n" % nYield)
    tex.write("Event yields at each step of the \\ttbar\\ selection %s. Figures in brackets do not account for the scale factors.\n" % tYield)
    tex.write("}\n")
    tex.write("\\begin{tabular}{l ||") #cut
    for process in procMC:
        tex.write("ll |") #each MC process
    tex.write("| ll || l}\n") #Total MC || Data
    tex.write("& ")
    for process in procMC:
        tex.write("\multicolumn{2}{c|")
        if process is procMC[len(procMC)-1]:
            tex.write("|")
        tex.write("}{%s} & " % process[0])
    tex.write("\\multicolumn{2}{c||}{Total MC} & Data \\\\ \n")
    tex.write("\\hline\\hline\n")
    hData = rootFile.Get("%s/%s" % (nYield, nYield))
    cut = 1
    while hData.GetBinContent(cut) > 0:
        tex.write("%s & " % (hData.GetXaxis().GetBinLabel(cut)).replace("#mu","$\\mu$").replace("+","$+$").replace(">","$>$"))
        wTot = 0
        uTot = 0
        for process in procMC:
            wMC =rootFile.Get("%s/%s_%s" % (nYield, nYield, process[1]))
            uMC =rootFile.Get("%s_no_weight/%s_no_weight_%s" % (nYield, nYield, process[1]))
            if abs(wMC.GetBinContent(cut) - uMC.GetBinContent(cut)) > 1e-6:
                tex.write(("$%.2e}$ & ($%.2e}$) & " % (wMC.GetBinContent(cut), uMC.GetBinContent(cut))).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))
            else :
                tex.write((" & ($%.2e}$) & " % uMC.GetBinContent(cut)).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))
            wTot = wTot + wMC.GetBinContent(cut) 
            uTot = uTot + uMC.GetBinContent(cut) 
        if abs(wTot - uTot) > 1e-6:
            tex.write(("$%.2e}$ & ($%.2e}$) & " % (wTot, uTot)).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))
        else :
            tex.write((" & ($%.2e}$) & " % uTot).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))
        tex.write(("$%.2e}$ \\\\ \n" % hData.GetBinContent(cut)).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-"))
        cut = cut+1
    tex.write("\\end{tabular}\n")
    tex.write("\\end{center}\n")
    tex.write("\\end{table}\n")

tex.write("\\clearpage\n")
tex.write("\n\\subsection{Efficiencies}\n")    

for hYield in hYields:
    tYield = hYield[0]
    nYield = hYield[1]
    tex.write("\\begin{table}[h!]\n")
    tex.write("\\begin{center}\n")
    tex.write("\\caption{\\label{%s}\n" % nYield)
    tex.write("Efficiency (\\%%) of each step of the \\ttbar\\ selection %s. Figures in brackets do not account for the scale factors.\n" % tYield)
    tex.write("}\n")
    tex.write("\\begin{tabular}{l ||") #cut
    for process in procMC:
        tex.write("rr |") #each MC process
    tex.write("| rr || r}\n") #Total MC || Data
    tex.write("& ")
    for process in procMC:
        tex.write("\multicolumn{2}{c|")
        if process is procMC[len(procMC)-1]:
            tex.write("|")
        tex.write("}{%s} & " % process[0])
    tex.write("\\multicolumn{2}{c||}{Total MC} & Data \\\\ \n")
    tex.write("\\hline\\hline\n")
    hData = rootFile.Get("%s/%s" % (nYield, nYield))
    cut = 2
    while hData.GetBinContent(cut) > 0:
        tex.write("%s & " % (hData.GetXaxis().GetBinLabel(cut)).replace("#mu","$\\mu$").replace("+","$+$").replace(">","$>$"))
        wNum = 0
        uNum = 0
        wDen = 0
        uDen = 0
        for process in procMC:
            wMC =rootFile.Get("%s/%s_%s" % (nYield, nYield, process[1]))
            uMC =rootFile.Get("%s_no_weight/%s_no_weight_%s" % (nYield, nYield, process[1]))
            if wMC.GetBinContent(cut-1) > 0 :
                if uMC.GetBinContent(cut-1) > 0 :
                    wEff = 100*wMC.GetBinContent(cut)/wMC.GetBinContent(cut-1)
                    uEff = 100*uMC.GetBinContent(cut)/uMC.GetBinContent(cut-1)
                    if abs(wEff - uEff) > 1e-15:
                        tex.write("$%.1f$ & ($%.1f$) & " % (wEff, uEff))
                    else :
                        tex.write(" & ($%.1f$) & " % uEff)
                else :        
                    wEff = 100*wMC.GetBinContent(cut)/wMC.GetBinContent(cut-1)
                    tex.write("$%.1f$ & -- & " % wEff)
            else :
                if uMC.GetBinContent(cut-1) > 0 :
                    uEff = 100*uMC.GetBinContent(cut)/uMC.GetBinContent(cut-1)
                    tex.write(" -- & $%.1f$ & " % uEff)
                else :        
                    tex.write(" -- & -- & ")
            wNum = wNum + wMC.GetBinContent(cut) 
            uNum = uNum + uMC.GetBinContent(cut) 
            wDen = wDen + wMC.GetBinContent(cut-1) 
            uDen = uDen + uMC.GetBinContent(cut-1) 
        if wDen > 0 :    
            if uDen > 0 :
                wTot = 100*wNum/wDen
                uTot = 100*uNum/uDen
                if abs(wTot - uTot) > 1e-15:
                    tex.write("$%.1f$ & ($%.1f$) & " % (wTot, uTot))
                else :
                    tex.write(" & ($%.1f$) & " % uTot)
            else :        
                wTot = 100*wNum/wDen
                tex.write("$%.1f$ & -- & " % wTot)
        else :        
            if uDen > 0 :
                uTot = 100*uNum/uDen
                tex.write(" -- & $%.1f$ & " % uTot)
            else :
                tex.write(" -- & -- & ")
        if hData.GetBinContent(cut-1) > 0 :    
            tex.write("$%.1f$ \\\\ \n" % (100*hData.GetBinContent(cut)/hData.GetBinContent(cut-1)))
        else :
            tex.write("-- \\\\ \n")
        cut = cut+1
    tex.write("\\end{tabular}\n")
    tex.write("\\end{center}\n")
    tex.write("\\end{table}\n")

tex.write("\\clearpage\n")
tex.write("\n\n\\section{Overall}")

for hDecay in hDecays:
    tDecay = hDecay[0]
    nDecay = hDecay[1]
    tex.write("\\begin{table}[h!]\n")
    tex.write("\\begin{center}\n")
    tex.write("\\caption{\\label{%s}\n" % nYield)
    tex.write("Event yields of a selection enriched in %s candidates. Figures in brackets do not account for the scale factors.\n" % tDecay)
    tex.write("}\n")
    tex.write("\\begin{tabular}{l ||") #cut
    hData = rootFile.Get("%s/%s" % (nDecay, nDecay))
    decay = 0
    while decay < hData.GetNbinsX():
        tex.write("ll |") #each decay channel
        decay = decay+1
        if decay == hData.GetNbinsX():
            tex.write("|")
    tex.write(" ll}\n") #all decay channels
    tex.write("& ")
    decay = 0
    while decay < hData.GetNbinsX():
        tex.write("\\multicolumn{2}{c|")
        if decay+1 == hData.GetNbinsX():
            tex.write("|")
        tex.write(("}{%s} & " % hData.GetXaxis().GetBinLabel(decay+1)).replace("#mu","$\\mu$").replace("+","$+$")) #each decay channel
        decay = decay+1
    tex.write("\\multicolumn{2}{c}{All} \\\\ \n") #all decay channels
    tex.write("\\hline\\hline\n")
    for process in procMC:
        tex.write("%s & " % process[0])
        wMC =rootFile.Get("%s/%s_%s" % (nDecay, nDecay, process[1]))
        uMC =rootFile.Get("%s_no_weight/%s_no_weight_%s" % (nDecay, nDecay, process[1]))
        decay = 1
        tW = 0
        tU = 0
        while decay < wMC.GetNbinsX()+1:
            tex.write(("$%.2e}$ & ($%.2e}$) & " % (wMC.GetBinContent(decay), uMC.GetBinContent(decay))).replace("0.00e+00","\\text{negligible}{").replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))
            tW = tW+wMC.GetBinContent(decay)
            tU = tU+uMC.GetBinContent(decay)
            decay = decay+1
        tex.write(("$%.2e}$ & ($%.2e}$) \\\\ \n " % (tW, tU)).replace("0.00e+00","\\text{negligible}{").replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))

    tex.write("\\hline\\hline\n")
    tex.write("Total MC & ")
    decay = 1
    tW = 0
    tU = 0
    while decay < hData.GetNbinsX()+1:
        wTot = 0
        uTot = 0
        for process in procMC:
            wMC =rootFile.Get("%s/%s_%s" % (nDecay, nDecay, process[1]))
            uMC =rootFile.Get("%s_no_weight/%s_no_weight_%s" % (nDecay, nDecay, process[1]))
            wTot = wTot + wMC.GetBinContent(decay)
            uTot = uTot + uMC.GetBinContent(decay)
        tex.write(("$%.2e}$ & ($%.2e}$) & " % (wTot, uTot)).replace("0.00e+00","\\text{negligible}{").replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))
        tW = tW+wTot
        tU = tU+uTot
        decay = decay+1
    tex.write(("$%.2e}$ & ($%.2e}$) \\\\ \n " % (tW, tU)).replace("0.00e+00","\\text{negligible}{").replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-").replace("e+","\\cdot10^{").replace("e-","\\cdot10^{-"))

    tex.write("\\hline\\hline\n")
    tex.write("Data & ")
    decay = 1
    tData = 0
    while decay < hData.GetNbinsX()+1:
        tex.write("\\multicolumn{2}{c|")
        if decay == hData.GetNbinsX():
            tex.write("|")
        tex.write(("}{$%.2e}$} & " % hData.GetBinContent(decay)).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-"))
        tData = tData+hData.GetBinContent(decay)
        decay = decay + 1
    tex.write(("\\multicolumn{2}{c}{$%.2e}$} \\\\ \n" % tData).replace("e+0","\\cdot10^{").replace("e-0","\\cdot10^{-"))

    tex.write("\\end{tabular}\n")
    tex.write("\\end{center}\n")
    tex.write("\\end{table}\n")


tex.write("\n%%%\n")
tex.write("\n\\end{document}")

tex.close()

os.chdir(texDir)
cmd = "pdflatex " + texFile
os.system(cmd)
cmd = "rm -f *.aux *.log"
os.system(cmd)
os.chdir("..")
print "\n"+os.path.join(texDir,texFile)+" has been created and compiled."
