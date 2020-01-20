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
import AtlasStyle, AtlasLabel
from collections import OrderedDict

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-m', '--mode', action='store', choices=['inclusive', 'per_category', 'per_process'], default='inclusive', help='run mode')
    parser.add_argument('--categorization', action='store', choices=['Event_XGB_16_Category', 'Event_XGB_QG_16_Category'], default='Event_XGB_16_Category', help='categorization')
    return parser.parse_args()

def Leg(hists, x0, y1):

    l=TLegend(x0, y1-0.05*len(hists), x0+0.46, y1)
    l.SetFillColor(0)
    l.SetLineColor(1)
    l.SetTextSize(0.04)
    l.SetShadowColor(0)
    l.SetFillStyle(0)
    l.SetBorderSize(0)
    for category in hists:
        l.AddEntry(hists[category], category)
    return l

def plot_var(samples, region, var, varname, nbin, minvar, maxvar, adjustNdivision=False, discrete=False, text_label=False, category=None, category_name=None, categorization='Event_XGB_16_Category'):

    basepath = 'skimmed_ntuples/'

    hist = OrderedDict()
    for i, sample_name in enumerate(samples):
        hist[sample_name] = TH1F('hist' + str(i), 'hist' + str(i), nbin, minvar, maxvar)

    for sample_name in samples:

        print sample_name

        for process in samples[sample_name]['process']:

            path = basepath + process

            for app in os.listdir(path):

                f = TFile(path + '/' + app)
                t = f.Get(region)
                htemp=TH1F('htemp','htemp', nbin, minvar, maxvar)
                htemp.Sumw2()
                t.Draw('%s>>htemp' % (var), '%s%s' %(samples[sample_name]['selection'], '' if not category else '*(%s==%d)' %(categorization, category)))
                hist[sample_name].Add(htemp)

    for sample_name in samples:
        hist[sample_name].Scale(1/hist[sample_name].Integral())
        hist[sample_name].SetLineColor(samples[sample_name]['color'])
        if 'width' in samples[sample_name]: hist[sample_name].SetLineWidth(samples[sample_name]['width'])
        if 'style' in samples[sample_name]: hist[sample_name].SetLineStyle(samples[sample_name]['style'])
        hist[sample_name].SetMarkerSize(0)

    C1 = TCanvas()
    C1.SetCanvasSize(400,370)
    C1.SetLeftMargin(0.15)
    C1.SetRightMargin(0.05)
    C1.SetTopMargin(0.05)
    C1.SetBottomMargin(0.15)

    THSvar = THStack(var,'')

    for sample_name in samples:
        THSvar.Add(hist[sample_name], 'hist E1')

    THSvar.Draw('NoStack')

    if adjustNdivision: THSvar.GetXaxis().SetNdivisions(adjustNdivision, kFALSE)

    THSvar.GetXaxis().SetTitle("%s" % (varname))
    #THSvar.GetXaxis().SetTitleOffset(1.5)
    THSvar.GetYaxis().SetTitle("Fraction of Events")
    THSvar.GetYaxis().SetTitleOffset(1.5)
    if discrete:
        for i in range(nbin):
            THSvar.GetXaxis().SetBinLabel(i+1, '%d' % minvar+i)
    if text_label:
         for i in range(nbin):
            THSvar.GetXaxis().SetBinLabel(i+1, text_label[i])
    max_value = THSvar.GetMaximum('NoStack')
    hist[samples.keys()[0]].GetYaxis().SetRangeUser(0, max_value * 1.4)

    AtlasLabel.ATLASLabel(0.185, 0.885, " Internal")
    AtlasLabel.myText(0.185, 0.885-0.063, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    AtlasLabel.myText(0.185,0.885-0.063*2, region.replace("zero_", "0-").replace("one_", "1-").replace("two_", "2-") if not category_name else category_name)

    l = Leg(hist, 0.52, 0.93)
    l.Draw("SAME")

    if not category:
        if not os.path.isdir('Variables/%s' % (region)):
            os.makedirs('Variables/%s' % (region))
        C1.Print("Variables/%s/%s.pdf" % (region, var.replace('/','_')))

    else:
        if not category_name: category_name = category

        if not os.path.isdir('%s/%s' % (categorization, category_name)):
            os.makedirs('%s/%s' % (categorization, category_name))
        C1.Print("%s/%s/%s.pdf" % (categorization, category_name, var.replace('/','_')))

    for sample_name in samples:
        hist[sample_name].Delete()

    return

def plotVars(region, category=None, category_name=None, categorization='Event_XGB_16_Category'):

    samples = OrderedDict([
        ('Data sideband', {'process': ['data_sid'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=180)', 'color': kBlack}),
        #('Data center', {'process': ['data_cen'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kBlue+3, 'style': 2}),
        ('VBF', {'process': ['VBF'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kRed, 'width': 3}),
        ('ggF', {'process': ['ggF'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kOrange+2, 'width': 3}),
        ('MC bkg sideband', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130))', 'color': kBlue}),
        ('MC bkg center', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kGreen+2})
    ])

    if not category:

        region_code = {'zero_jet': 0, 'one_jet': 1, 'two_jet': 2}

        plot_var(samples, region,'ClassOut_XGB_Higgs', 'O_{ggF}^{(%d)}' % region_code[region], 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
        if region == 'two_jet': plot_var(samples, region,'ClassOut_XGB_VBF', 'O_{VBF}', 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region,'ClassOut_XGB_QG_Higgs', 'O_{ggF}^{(%d)}' % region_code[region], 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
        if region == 'two_jet': plot_var(samples, region,'ClassOut_XGB_QG_VBF', 'O_{VBF}', 20, 0, 1, category=category, category_name=category_name, categorization=categorization)

    samples = OrderedDict([
        ('Data sideband', {'process': ['data_sid'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=180)', 'color': kBlack}),
        #('Data center', {'process': ['data_cen'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kBlue+3, 'style': 2}),
        ('Higgs signal', {'process': ['VBF', 'VH', 'ggF', 'ttH'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kRed, 'width': 3}),
        ('MC bkg sideband', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130))', 'color': kBlue}),
        ('MC bkg center', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': kGreen+2})
    ])

    plot_var(samples, region,'m_mumu', 'm_{#mu#mu} [GeV]', 20, 110, 180, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'Z_PT_OnlyNearFsr', 'p_{T}^{#mu#mu} [GeV]', 20, 0, 200 if 'two_jet' in region else 100, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'Z_Y_OnlyNearFsr', 'y_{#mu#mu}', 20, -4, 4, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'fabs(Muons_CosThetaStar)', '|cos#theta*|', 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'Muons_CosThetaStar', 'cos#theta*', 20, -1, 1, category=category, category_name=category_name, categorization=categorization)

    if region != 'zero_jet':
        plot_var(samples, region, 'Jets_PT_Lead', 'p_{T}^{j_{#scale[1]{1}}} [GeV]', 20, 0, 300, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Eta_Lead', '#eta_{j_{1}}', 20, -5, 5, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'DeltaPhi_mumuj1', '#Delta#phi_{#mu#mu,j_{1}}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_QGscore_Lead', 'N_{tracks}^{j_{1}}', 20, 0, 20, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_QGflag_Lead', 'QG_{j_{1}}', 3, -1, 2, text_label=['Untagged', 'Quark', 'Gluon'], category=category, category_name=category_name, categorization=categorization)

    if region == 'two_jet':
        plot_var(samples, region, 'Jets_PT_Sub', 'p_{T}^{j_{#scale[1]{2}}} [GeV]', 20, 0, 200, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Eta_Sub', '#eta_{j_{2}}', 20, -4.5, 4.5, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'DeltaPhi_mumuj2', '#Delta#phi_{#mu#mu,j_{2}}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_QGscore_Sub', 'N_{tracks}^{j_{2}}', 20, 0, 20, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_QGflag_Sub', 'QG_{j_{2}}', 3, -1, 2, text_label=['Untagged', 'Quark', 'Gluon'], category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_PT_jj', 'p_{T}^{jj} [GeV]', 20, 0, 400, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Y_jj', 'y_{jj}', 20, -6, 6, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'DeltaPhi_mumujj', '#Delta#phi_{#mu#mu,jj}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Minv_jj', 'm_{jj} [GeV]' , 20, 0, 1000, adjustNdivision = 40205, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region,'Event_MET', 'E_{T}^{miss} [GeV]', 20, 0, 150 if 'two_jet' in region else 80, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region,'DeltaPhi_mumuMET', '#Delta#phi_{#mu#mu,MET}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region,'Event_Ht-Muons_PT_Lead-Muons_PT_Sub', 'H_{T} [GeV]', 20, 0, 200, category=category, category_name=category_name, categorization=categorization)

def main():

    args=getArgs()

    gROOT.SetBatch(True)
    TH1.SetDefaultSumw2(1)

    if args.mode == 'inclusive':

        regions = ['zero_jet', 'one_jet', 'two_jet']

        for region in regions:
        
            plotVars(region)

    elif args.mode == 'per_category':

        categories = {
            1: 'VBF-High',
            2: 'VBF-Medium',
            3: 'VBF-Low',
            4: 'VBF-VeryLow',
            5: 'Higgs-2Jet-High',
            6: 'Higgs-2Jet-Medium',
            7: 'Higgs-2Jet-Low',
            8: 'Higgs-2Jet-VeryLow',
            9: 'Higgs-1Jet-High',
            10: 'Higgs-1Jet-Medium',
            11: 'Higgs-1Jet-Low',
            12: 'Higgs-1Jet-VeryLow',
            13: 'Higgs-0Jet-High',
            14: 'Higgs-0Jet-Medium',
            15: 'Higgs-0Jet-Low',
            16: 'Higgs-0Jet-VeryLow',
        }

        for category in categories:

            if category <= 8: region = 'two_jet'
            elif category <= 12: region = 'one_jet'
            else: region = 'zero_jet'

            print "Process category: ", categories[category]

            plotVars(region, category=category, category_name=categories[category], categorization=args.categorization)

    return

if __name__ == '__main__':
    main()
