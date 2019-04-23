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
from ROOT import *
from argparse import ArgumentParser
from AtlasStyle import *

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['zero_jet', 'one_jet', 'two_jet', 'all_jet'], default='zero_jet', help='Region to process')
    return  parser.parse_args()

def SetStyle(hist,i):

    FillColor=[kGreen-9,kBlue-9,kBlue+3,kRed-9,0,0,0]
    LineColor=[kGreen+3,kBlue,kBlue+3,kRed+2,kOrange,kGreen,1]
    hist.SetFillColorAlpha(FillColor[i],0.25)
    hist.SetLineColor(LineColor[i])
    hist.SetMarkerSize(0)
    #if i==2: hist.SetFillStyle(3554)
    
    return

def Leg(hists,x0,y1):

    n=len(hists)
    l=TLegend(x0,y1-n*0.06,x0+0.46,y1)
    l.SetFillColor(0)
    l.SetLineColor(1)
    l.SetTextSize(0.04)
    l.SetShadowColor(0)
    l.SetFillStyle(0)
    l.SetBorderSize(0)
    l.AddEntry(hists['bkg'],'Data sideband')
    l.AddEntry(hists['ggF'],'ggF')
    l.AddEntry(hists['VBF'],'VBF')
#    l.AddEntry(hists['sig'],'ggF + VBF')
    return l

def plot_var(categories,region,train_sig,branch,var,lend,rend):

    hist={}
    exist={}
    i=0
    for category in categories:
        hist[category]=TH1F('hist_%s_%s'%(var,category),'hist_%s_%s'%(var,category),40,lend,rend)

    for category in categories:
        exist[category]=False
        for app in os.listdir('outputs/model_%s%s'%(train_sig, region)):
            
            if category not in app: continue
            if ".root" not in app: continue

            #print app
            exist[category]=True

            f=TFile('outputs/model_%s%s/%s'%(train_sig,region,app))
            t=f.Get('test')
            
            htemp=TH1F('htemp','htemp',40,lend,rend)

            #if 'mc16a' in app: lumi=36
            #elif 'mc16d' in app: lumi=44
            t.Draw("%s>>htemp"%(branch),"weight")


            hist[category]+=htemp


        i+=1

    #print hist
    #hist['sig'].Scale(1/hist['sig'].Integral())
    #hist_bkg_sid.SetFillColorAlpha(kBlue+3, 0.25)
#    hist['sig'].SetLineColor(kRed)
#    hist['sig'].SetMarkerSize(0)
    #hist_bkg_sid.SetFillStyle(3554)

    hist['ggF'].SetLineColor(kGreen+3)
    hist['ggF'].SetMarkerSize(0)

    hist['VBF'].SetLineColor(kRed)
    hist['VBF'].SetMarkerSize(0)


    #hist['bkg'].Scale(1/hist['bkg'].Integral())
    #hist_bkg_sid.SetFillColorAlpha(kBlue+3, 0.25)
    hist['bkg'].SetLineColor(kBlue)
    hist['bkg'].SetMarkerSize(0)
    #hist_bkg_sid.SetFillStyle(3554)

    C1=TCanvas()
    #C1.Range(0,0,1,1)
    C1.SetCanvasSize(400,370)
    C1.SetLeftMargin(0.135)
    C1.SetRightMargin(0.05)
    C1.SetTopMargin(0.05)
    C1.SetBottomMargin(0.135)

    THSvar=THStack(var,'')

    for category in categories:
        if not exist[category]: continue
        if hist[category].Integral() != 0:
            hist[category].Scale(1/hist[category].Integral())

            THSvar.Add(hist[category],"HIST E1")

    THSvar.Draw("Nostack")

    max_value = THSvar.GetMaximum('NoStack')
    hist['bkg'].GetYaxis().SetRangeUser(0, max_value * 1.4)

    THSvar.GetXaxis().SetTitle(var)
    #jet1_pt.GetXaxis().SetTitleOffset(1.5)
    THSvar.GetYaxis().SetTitle("1/N dN/dX")
    #jet1_pt.GetYaxis().SetTitleOffset(1.5)
    AtlasLabel.ATLASLabel(0.17,0.88," Internal")
    AtlasLabel.myText(0.17,0.88-0.063,"#sqrt{s} = 13 TeV, 140 fb^{-1}")

    l=Leg(hist,0.56,0.92)
    l.Draw("SAME")


    if not os.path.isdir('bdt_plots'):
        os.makedirs('bdt_plots')
    C1.Print('bdt_plots/model_%s%s%s%s.pdf'%(train_sig,region, '_VBF' if '_VBF' in branch else '', '_t' if '_t' in branch else ''))

    return


def main():

    gROOT.SetBatch(True)

    args=getArgs()

    region = args.region

    bkgs=['bkg']
    sigs=['sig']
    vbf=['VBF']
    ggF=['ggF']

    train_sig=''

    #plotVars(bkgs+sigs,region,train_sig)
    #plot_var(bkgs+sigs,region,train_sig,'bdt_score_t','O_{BDT}',0,1)
    plot_var(bkgs+vbf+ggF,region,train_sig,'bdt_score','O_{BDT}',0,1)
    plot_var(bkgs+vbf+ggF,region,train_sig,'bdt_score_t','O_{BDT}',0,1)
    if region == 'two_jet' or region == 'all_jet': 
        plot_var(bkgs+vbf+ggF,region,train_sig,'bdt_score_VBF','O_{BDT}',0,1)
        plot_var(bkgs+vbf+ggF,region,train_sig,'bdt_score_VBF_t','O_{BDT}',0,1)

    return

if __name__ == '__main__':
    main()
