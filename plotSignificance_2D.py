#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.28.2018
#
#
#
#
#
import os
from ROOT import TH2F, TCanvas, gStyle, gROOT
from argparse import ArgumentParser
from AtlasStyle import *
from math import sqrt

def main():

    gROOT.SetBatch(True)

    min_nggf = 2
    max_nggf = 7
    min_nvbf = 2
    max_nvbf = 7

    sigplot=TH2F('impplot','impplot', max_nggf - min_nggf + 1, min_nggf, max_nggf + 1, max_nvbf - min_nvbf + 1, min_nvbf, max_nvbf + 1)

    # test set
    #zs=[
    #    [2, 7, 1.381], [3, 7, 1.401], [4, 7, 1.397], [5,  7, 1.412], [6,  7, 1.420], [7,  7, 1.423]#, [8,  7, 1.458], [9,  7, 1.451]
    #   ,[2, 6, 1.390], [3, 6, 1.410], [4, 6, 1.409], [5,  6, 1.410], [6,  6, 1.415], [7,  6, 1.418]#, [8,  6, 1.448], [9,  6, 1.438]
    #   ,[2, 5, 1.379], [3, 5, 1.389], [4, 5, 1.409], [5,  5, 1.419], [6,  5, 1.424], [7,  5, 1.427]#, [8,  5, 1.448], [9,  5, 1.438]
    #   ,[2, 4, 1.358], [3, 4, 1.378], [4, 4, 1.418], [5,  4, 1.408], [6,  4, 1.410], [7,  4, 1.404]#, [8,  4, 1.449], [9,  4, 1.453]
    #   ,[2, 3, 1.379], [3, 3, 1.405], [4, 3, 1.382], [5,  3, 1.389], [6,  3, 1.401], [7,  3, 1.402]#, [8,  3, 1.440], [9,  3, 1.448]
    #   ,[2, 2, 1.331], [3, 2, 1.383], [4, 2, 1.406], [5,  2, 1.412], [6,  2, 1.418], [7,  2, 1.420]#, [8,  2, 1.408], [9,  2, 1.408]
    #   ]

    # cate set
    zs=[
        [2, 7, 1.404], [3, 7, 1.424], [4, 7, 1.433], [5,  7, 1.439], [6,  7, 1.442], [7,  7, 1.444]#, [8,  7, 1.458], [9,  7, 1.451]
       ,[2, 6, 1.402], [3, 6, 1.422], [4, 6, 1.432], [5,  6, 1.437], [6,  6, 1.440], [7,  6, 1.443]#, [8,  6, 1.448], [9,  6, 1.438]
       ,[2, 5, 1.398], [3, 5, 1.419], [4, 5, 1.429], [5,  5, 1.435], [6,  5, 1.438], [7,  5, 1.441]#, [8,  5, 1.448], [9,  5, 1.438]
       ,[2, 4, 1.392], [3, 4, 1.414], [4, 4, 1.425], [5,  4, 1.431], [6,  4, 1.435], [7,  4, 1.437]#, [8,  4, 1.449], [9,  4, 1.453]
       ,[2, 3, 1.382], [3, 3, 1.407], [4, 3, 1.418], [5,  3, 1.426], [6,  3, 1.430], [7,  3, 1.432]#, [8,  3, 1.440], [9,  3, 1.448]
       ,[2, 2, 1.355], [3, 2, 1.390], [4, 2, 1.408], [5,  2, 1.417], [6,  2, 1.422], [7,  2, 1.425]#, [8,  2, 1.408], [9,  2, 1.408]
       ]

    for z in zs:
        sigplot.Fill(z[0], z[1], z[2])

    if not os.path.isdir('significances/sigplots'):
        print 'INFO: Creating new folder: \"significances/sigplots\"'
        os.makedirs("significances/sigplots")


    C2=TCanvas()
    #C2.Range(0,0,1,1)
    C2.SetCanvasSize(550,395)
    C2.SetLeftMargin(0.115)
    C2.SetRightMargin(0.194)
    C2.SetTopMargin(0.15)
    C2.SetBottomMargin(0.135)
    gStyle.SetPaintTextFormat("4.3f")

    sigplot.Draw("COLZTEXT")
    sigplot.GetXaxis().SetTitle("Number of ggF bins")
    sigplot.GetXaxis().SetTitleOffset(1.2)
    sigplot.GetYaxis().SetTitle("Number of VBF bins")
    sigplot.GetYaxis().SetTitleOffset(0.9)
    sigplot.GetZaxis().SetTitle("Significance [#sigma]")
    sigplot.GetZaxis().SetTitleOffset(1.4)
    sigplot.SetMarkerSize(2.5)
    sigplot.SetAxisRange(1.330, 1.444,"Z")

    xlabel=range(min_nggf, max_nggf + 1)
    ylabel=range(min_nvbf, max_nvbf + 1)

    for i in range(max_nggf - min_nggf + 1):
        sigplot.GetXaxis().SetBinLabel(i+1,'%d' % xlabel[i])
    for i in range(max_nvbf - min_nvbf + 1):
        sigplot.GetYaxis().SetBinLabel(i+1,'%d' % ylabel[i])

    AtlasLabel.ATLASLabel(0.135,0.95," Internal")
    AtlasLabel.myText(0.37,0.95,"#sqrt{s} = 13 TeV, 140 fb^{-1}")
    AtlasLabel.myText(0.135,0.88,"120 GeV < m_{#mu#mu} < 130 GeV")
    #C2.Print("significances/sigplots/ForCONF_2D.pdf")
    C2.Print("significances/sigplots/ForCONF_2D_catset.pdf")

    return

if __name__ == '__main__':
    main()
