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

def Leg(hists,x0,y1,region):

    l=TLegend(x0,y1-0.05*len(hists),x0+0.46,y1)
    l.SetFillColor(0)
    l.SetLineColor(1)
    l.SetTextSize(0.04)
    l.SetShadowColor(0)
    l.SetFillStyle(0)
    l.SetBorderSize(0)

    for i in range(len(hists)):
        l.AddEntry(hists[i], region+' %d'%(i+1))
    return l

def plot_var(category, regions, region,var, varname, nbin,minvar,maxvar):


    f = TFile('outputs_2D/model_%s/%s.root'%(region if region != 'two_jet_VBF' else 'two_jet', category))
    t = f.Get('test')

    colors = [kRed, kGreen+2, kBlue, kViolet, kBlack]

    hist = []
    for i in range(1, regions[region] + 1):
        hist.append(TH1F('hist_%d'%(i), 'hist_%d'%(i), nbin, minvar,maxvar))
        hist[i-1].Sumw2()
        t.Draw('%s>>hist_%d' % (var, i),'weight*(m_mumu>=110&&m_mumu<=180&&bdt_category==%d)'%(i if region != 'two_jet' else i + regions['two_jet_VBF']))
        if hist[i-1].Integral() != 0: hist[i-1].Scale(1/hist[i-1].Integral())
        else: print 'Integral = 0!!!!! ', i
        hist[i-1].SetLineColor(colors[i-1])
        hist[i-1].SetMarkerSize(0.5)
        hist[i-1].SetMarkerColor(colors[i-1])


    C1=TCanvas()
    #C1.Range(0,0,1,1)
    C1.SetCanvasSize(400,370)
    C1.SetLeftMargin(0.135)
    C1.SetRightMargin(0.05) #0.05
    C1.SetTopMargin(0.05)
    C1.SetBottomMargin(0.142)

    THSvar=THStack(var,'')

    for i in range(1, regions[region] + 1):
        THSvar.Add(hist[i-1],'hist E1')

    THSvar.Draw('NoStack')

    #THSvar.GetXaxis().SetNdivisions(40205, kFALSE)

    THSvar.GetXaxis().SetTitle("%s"%(varname))
    #jet1_pt.GetXaxis().SetTitleOffset(1.5)
    THSvar.GetYaxis().SetTitle("1/N dN/dX")
    #jet1_pt.GetYaxis().SetTitleOffset(1.5)
    #if var[0] == 'n' or var[0] == 'N':
    #    for i in range(3 if region == 'two_jet' else 2):
    #        THSvar.GetXaxis().SetBinLabel(i+1,'%d' % i)
    max_value = THSvar.GetMaximum('NoStack')
    hist[0].GetYaxis().SetRangeUser(0, max_value * 1.4)

    AtlasLabel.ATLASLabel(0.17,0.88," Internal")
    AtlasLabel.myText(0.17,0.88-0.063,"#sqrt{s} = 13 TeV, 140 fb^{-1}")

    l=Leg(hist,0.56,0.93,region)
    l.Draw("SAME")


    if not os.path.isdir('plots/%s'%(region)):
        os.makedirs('plots/%s'%(region))
    C1.Print("plots/%s/%s_%s.pdf"%(region,category,var.replace('/','_')))

    
    

    return


def main():

    gROOT.SetBatch(True)
    #TGaxis.SetMaxDigits(1)

    sigs = ['VBF', 'VH', 'ggF', 'ttH']
    bkgs = ['data']
    #regions = ['zero_jet', 'one_jet', 'two_jet', 'two_jet_VBF']
    regions = {'zero_jet': 5, 'one_jet': 5, 'two_jet': 5, 'two_jet_VBF': 5}


    for category in sigs+bkgs+['sig']:

        print category

        for region in regions:

            print region
            plot_var(category, regions, region,'Muons_Minv_MuMu_Sigma/m_mumu', 'Fit Sigma/m_{#mu#mu}',20,0,0.04 if 'two_jet' in region else 0.03)
            plot_var(category, regions, region,'Muons_MCPExpReso_Smear_Minv_MuMu_Sigma/m_mumu', 'MCP Sigma/m_{#mu#mu}',20,0.01,0.025 if region == 'zero_jet' else 0.03)
            #continue

            plot_var(category, regions, region,'m_mumu', 'm_{#mu#mu} [GeV]', 20, 110, 180)

            plot_var(category, regions, region,'metFinalTrk', 'p_{T}^{miss} [GeV]', 20, 0, 150 if 'two_jet' in region else 80)
            #plot_var(sigs, bkgs, region,'metFinalTrkPhi',20,-3.14,3.14)
            plot_var(category, regions, region,'Z_PT', 'p_{T#mu#mu} [GeV]', 20,0,200 if 'two_jet' in region else 100)
            plot_var(category, regions, region,'Z_Eta', '#eta_{#mu#mu} [GeV]', 20,-4,4)
            #plot_var(category, regions, region,'Z_Phi', '#phi_{#mu#mu}', 20,-3.14,3.14)
            plot_var(category, regions, region,'fabs(Muons_CosThetaStar)', '|cos#theta^{*}|', 20,0,1)
            plot_var(category, regions, region,'Muons_Minv_MuMu_Sigma', 'Fit Sigma [GeV]',20,0,4 if region == 'zero_jet' else 5)
            plot_var(category, regions, region,'Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', 'MCP Sigma [GeV]',20,1,5)
            if region == 'two_jet' or region == 'two_jet_VBF':
                plot_var(category, regions, region,'Jets_Minv_jj', 'm_{jj} [GeV]' ,20,0,1000)
                plot_var(category, regions, region,'Jets_PT_jj', 'p_{Tjj} [GeV]', 20,0,400)
                plot_var(category, regions, region,'Jets_Eta_jj', '#eta_{jj}', 20,-6,6)
                plot_var(category, regions, region,'DeltaPhi_mumujj', '#Delta#phi_{#mu#mu,jj}', 20,-3.14,3.14)
            if region != 'zero_jet':
                plot_var(category, regions, region,'Jets_PT_Lead', 'p_{Tj1} [GeV]', 20,0,300)
                plot_var(category, regions, region,'Jets_Eta_Lead', '#eta_{j1}', 20,-5,5)
                plot_var(category, regions, region,'DeltaPhi_mumuj1', '#Delta#phi_{#mu#mu,j1}', 20,-3.14,3.14)
            if region == 'two_jet' or region == 'two_jet_VBF':
                plot_var(category, regions, region,'Jets_PT_Sub', 'p_{Tj2} [GeV]', 20,0,200)
                plot_var(category, regions, region,'Jets_Eta_Sub', '#eta_{j2}', 20,-4.5,4.5)
                plot_var(category, regions, region,'DeltaPhi_mumuj2', '#Delta#phi_{#mu#mu,j2}', 20,-3.14,3.14)

    return

if __name__ == '__main__':
    main()
