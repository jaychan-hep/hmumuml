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
import ROOT
from argparse import ArgumentParser
import AtlasStyle, AtlasLabel
from collections import OrderedDict

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-m', '--mode', action='store', choices=['inclusive', 'per_category', 'per_process', "m_mumu", "BDT_sideband_center_ratio"], default='inclusive', help='run mode')
    parser.add_argument('--categorization', action='store', choices=['Event_Paper_Category'], default='Event_Paper_Category', help='categorization')
    return parser.parse_args()

def Leg(hists, x0, y1):

    l = ROOT.TLegend(x0, y1-0.05*len(hists), x0+0.46, y1)
    l.SetFillColor(0)
    l.SetLineColor(1)
    l.SetTextSize(0.035)
    l.SetShadowColor(0)
    l.SetFillStyle(0)
    l.SetBorderSize(0)
    for category in hists:
        l.AddEntry(hists[category], category)
    return l

def plot_var_per_sample(samples, region, categories, var, varname, nbin, minvar, maxvar, adjustNdivision=False, discrete=False, logY=False, text_label=False, sample_name=None, region_name=None, categorization='Event_Paper_Category', save_histogram=False):

    print('==============================================')
    print('INFO: Samples ->', samples)
    print('INFO: Variable ->', var)

    basepath = 'skimmed_ntuples/'

    hist = OrderedDict()
    for i, category_name in enumerate(categories):
        hist[category_name] = ROOT.TH1F(('hist' + str(i)) if not 'save string' in categories[category_name] else categories[category_name]['save string'], ('hist' + str(i)) if not 'save string' in categories[category_name] else categories[category_name]['save string'], nbin, minvar, maxvar)

    for category_name in categories:

        print('INFO: Include category ', category_name)

        for sample in samples:

            path = basepath + sample

            for app in os.listdir(path):

                f = ROOT.TFile(path + '/' + app)
                t = f.Get(region)
                htemp = ROOT.TH1F('htemp','htemp', nbin, minvar, maxvar)
                htemp.Sumw2()
                t.Draw('%s>>htemp' % (var), '%s%s' %(categories[category_name]['selection'], '' if 'category' not in categories[category_name] else '*(%s==%d)' %(categorization, categories[category_name]['category'])))
                hist[category_name].Add(htemp)

    for category_name in categories:
        if hist[category_name].Integral() != 0: hist[category_name].Scale(1/hist[category_name].Integral())
        hist[category_name].SetLineColor(categories[category_name]['color'])
        if 'width' in categories[category_name]: hist[category_name].SetLineWidth(categories[category_name]['width'])
        if 'style' in categories[category_name]: hist[category_name].SetLineStyle(categories[category_name]['style'])
        hist[category_name].SetMarkerSize(0)

    C1 = ROOT.TCanvas()
    C1.SetCanvasSize(400,370)
    C1.SetLeftMargin(0.15)
    C1.SetRightMargin(0.05)
    C1.SetTopMargin(0.05)
    C1.SetBottomMargin(0.15)
    if logY: C1.SetLogy()

    THSvar = ROOT.THStack(var,'')

    for category_name in categories:
        THSvar.Add(hist[category_name], 'hist E1')

    THSvar.Draw('NoStack')

    if adjustNdivision: THSvar.GetXaxis().SetNdivisions(adjustNdivision, ROOT.kFALSE)

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
    if logY:
        min_value = THSvar.GetMinimum('NoStack')
        if max_value/min_value > 10**3: min_value = max_value/10**3
        hist[list(categories)[0]].GetYaxis().SetRangeUser(min_value, max_value**1.4/min_value**0.4)
    else:
        hist[list(categories)[0]].GetYaxis().SetRangeUser(0, max_value * 1.4)

    AtlasLabel.ATLASLabel(0.185, 0.885, " Internal")
    AtlasLabel.myText(0.185, 0.885-0.063, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    if sample_name: AtlasLabel.myText(0.185,0.885-0.063*2, sample_name)
    else: sample_name = "sample"

    l = Leg(hist, 0.52, 0.93)
    l.Draw("SAME")

    if not region_name:
        region_name = region

    if not os.path.isdir('%s/%s' % (categorization, sample_name)):
        os.makedirs('%s/%s' % (categorization, sample_name))
    C1.Print("%s/%s/%s.pdf" % (categorization, sample_name, var.replace('/','_') + '_' + region_name))

    if save_histogram:
        outfile = ROOT.TFile("%s/%s/%s.root" % (categorization, sample_name, var.replace('/','_') + '_' + region_name), "RECREATE")
        for category_name in categories:
            hist[category_name].Write()

    for category_name in categories:
        hist[category_name].Delete()

    return

def plot_var(samples, region, var, varname, nbin, minvar, maxvar, adjustNdivision=False, discrete=False, logY=False, text_label=False, category=None, category_name=None, categorization='Event_Paper_Category', region_mark=None, arrow=None):

    print('==============================================')
    if category:
        if not category_name: category_name = category
        print('INFO: Processing region ', category_name)
    else:
        print('INFO: Processing region ', region)
    print('INFO: Variable ->', var)

    basepath = 'skimmed_ntuples/'

    hist = OrderedDict()
    for i, sample_name in enumerate(samples):
        hist[sample_name] = ROOT.TH1F('hist' + str(i), 'hist' + str(i), nbin, minvar, maxvar)

    for sample_name in samples:

        print('INFO: Include sample ', sample_name)

        for process in samples[sample_name]['process']:

            path = basepath + process

            for app in os.listdir(path):

                f = ROOT.TFile(path + '/' + app)
                t = f.Get(region)
                htemp = ROOT.TH1F('htemp','htemp', nbin, minvar, maxvar)
                htemp.Sumw2()
                t.Draw('%s>>htemp' % (var), '%s%s' %(samples[sample_name]['selection'], '' if not category else '*(%s==%d)' %(categorization, category)))
                hist[sample_name].Add(htemp)

    for sample_name in samples:
        if hist[sample_name].Integral() != 0: hist[sample_name].Scale(1/hist[sample_name].Integral())
        hist[sample_name].SetLineColor(samples[sample_name]['color'])
        if 'width' in samples[sample_name]: hist[sample_name].SetLineWidth(samples[sample_name]['width'])
        if 'style' in samples[sample_name]: hist[sample_name].SetLineStyle(samples[sample_name]['style'])
        hist[sample_name].SetMarkerSize(0)

    C1 = ROOT.TCanvas()
    C1.SetCanvasSize(400,370)
    C1.SetLeftMargin(0.15)
    C1.SetRightMargin(0.05)
    C1.SetTopMargin(0.05)
    C1.SetBottomMargin(0.15)
    if logY: C1.SetLogy()

    THSvar = ROOT.THStack(var,'')

    for sample_name in samples:
        THSvar.Add(hist[sample_name], 'hist E1')

    THSvar.Draw('NoStack')

    if adjustNdivision: THSvar.GetXaxis().SetNdivisions(adjustNdivision, ROOT.kFALSE)

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
    if logY:
        min_value = THSvar.GetMinimum('NoStack')
        if max_value/min_value > 10**3: min_value = max_value/10**3
        hist[list(samples)[0]].GetYaxis().SetRangeUser(min_value, max_value**1.6/min_value**0.6)
    else:
        hist[list(samples)[0]].GetYaxis().SetRangeUser(0, max_value * 1.4)

    AtlasLabel.ATLASLabel(0.185, 0.885, "")
    AtlasLabel.myText(0.185, 0.885-0.063, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    AtlasLabel.myText(0.185,0.885-0.063*2, region.replace("zero_", "0-").replace("one_", "1-").replace("two_", "2-") if not category_name else category_name)

    l = Leg(hist, 0.52, 0.93)
    l.Draw("SAME")

    if region_mark:
        region_mark_line = {}
        for region_boundary in region_mark:
            region_mark_line[region_boundary] = ROOT.TLine(region_boundary, 0, region_boundary, max_value * 0.6 if not logY else max_value**0.9/min_value**(-0.1))
            region_mark_line[region_boundary].SetLineStyle(2)
            region_mark_line[region_boundary].SetLineColor(ROOT.kGray + 1)
            region_mark_line[region_boundary].SetLineWidth(2)
            region_mark_line[region_boundary].Draw()

    if arrow:
        arrow_with_txt = {}
        for ar in arrow:
            arrow_with_txt[ar] = ROOT.TArrow(arrow[ar][0], max_value * 0.55 if not logY else max_value * 0.3, arrow[ar][1], max_value * 0.55 if not logY else max_value * 0.3, 0.05, "|>")
            arrow_with_txt[ar].SetAngle(40)
            arrow_with_txt[ar].SetLineWidth(2)
            arrow_with_txt[ar].SetLineColor(ROOT.kGray + 1)
            arrow_with_txt[ar].SetFillColor(ROOT.kGray + 1)
            arrow_with_txt[ar].Draw()

            arrow_txt = ROOT.TLatex(0.7*arrow[ar][0] + 0.3*arrow[ar][1], 0.5 if not logY else 0.65, ar)
            arrow_txt.SetNDC()
            # arrow_txt.SetTextAngle()
            arrow_txt.SetTextSize(0.04)
            arrow_txt.SetTextColor(ROOT.kGray + 1)
            arrow_txt.Draw()

    if not category:
        if not os.path.isdir('Variables/%s' % (region)):
            os.makedirs('Variables/%s' % (region))
        C1.Print("Variables/%s/%s.pdf" % (region, var.replace('/','_')))

    else:

        if not os.path.isdir('%s/%s' % (categorization, category_name)):
            os.makedirs('%s/%s' % (categorization, category_name))
        C1.Print("%s/%s/%s.pdf" % (categorization, category_name, var.replace('/','_')))

    for sample_name in samples:
        hist[sample_name].Delete()

    return

def plot_var_ratio(samples, region, var, varname, nbin, minvar, maxvar, adjustNdivision=False, discrete=False, text_label=False, category=None, category_name=None, ratio_title=None, categorization='Event_Paper_Category'):

    print('==============================================')
    if category:
        if not category_name: category_name = category
        print('INFO: Processing region ', category_name)
    else:
        print('INFO: Processing region ', region)
    print('INFO: Variable ->', var)

    basepath = 'skimmed_ntuples/'

    hist = OrderedDict()
    for i, sample_name in enumerate(samples):
        hist[sample_name] = ROOT.TH1F('hist' + str(i), 'hist' + str(i), nbin, minvar, maxvar)

    hist_ratio = OrderedDict()
    denominator = OrderedDict()
    for sample in samples:
        if 'ratio' not in samples[sample]: continue
        if samples[sample]['ratio'] == 'denominator':
            denominator['all'] = sample
        elif samples[sample]['ratio'] == 'numerator':
            hist_ratio[sample] = ''
        else:
            denominator[samples[sample]['ratio']] = sample

    for sample_name in samples:

        print('INFO: Include sample ', sample_name)

        for process in samples[sample_name]['process']:

            path = basepath + process

            for app in os.listdir(path):

                f = ROOT.TFile(path + '/' + app)
                t = f.Get(region)
                htemp = ROOT.TH1F('htemp','htemp', nbin, minvar, maxvar)
                htemp.Sumw2()
                t.Draw('%s>>htemp' % (var), '%s%s' %(samples[sample_name]['selection'], '' if not category else '*(%s==%d)' %(categorization, category)))
                hist[sample_name].Add(htemp)

    for sample_name in samples:
        if hist[sample_name].Integral() != 0: hist[sample_name].Scale(1/hist[sample_name].Integral())
        hist[sample_name].SetLineColor(samples[sample_name]['color'])
        if 'width' in samples[sample_name]: hist[sample_name].SetLineWidth(samples[sample_name]['width'])
        if 'style' in samples[sample_name]: hist[sample_name].SetLineStyle(samples[sample_name]['style'])
        hist[sample_name].SetMarkerSize(0)

    for sample in hist_ratio:
        hist_ratio[sample] = hist[sample].Clone()
        htemp = hist[denominator['all'] if 'all' in denominator else denominator[sample]].Clone()
        hist_ratio[sample].Divide(htemp)

    C1 = ROOT.TCanvas()
    C1.SetCanvasSize(800,800)

    pad1 = ROOT.TPad("pad1", "pad2", 0, 0.3, 1, 1.0)
    #pad1.SetGridx()
    pad1.SetTopMargin(0.05)
    pad1.SetLeftMargin(0.15)
    pad1.Draw()
    pad1.cd()

    THSvar = ROOT.THStack(var,'')

    for sample_name in samples:
        THSvar.Add(hist[sample_name], 'hist E1')

    THSvar.Draw('NoStack')

    if adjustNdivision: THSvar.GetXaxis().SetNdivisions(adjustNdivision, ROOT.kFALSE)

    THSvar.GetXaxis().SetLabelSize(0)
    THSvar.GetYaxis().SetLabelSize(26)
    THSvar.GetYaxis().SetLabelFont(43)

    axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axis.SetLabelFont(43)
    axis.SetLabelSize(15)
    axis.Draw()

    THSvar.GetYaxis().SetTitle("Fraction of Events")
    THSvar.GetYaxis().SetTitleSize(26)
    THSvar.GetYaxis().SetTitleFont(43)
    THSvar.GetYaxis().SetTitleOffset(1.8)

    if discrete:
        for i in range(nbin):
            THSvar.GetXaxis().SetBinLabel(i+1, '%d' % minvar+i)
    if text_label:
         for i in range(nbin):
            THSvar.GetXaxis().SetBinLabel(i+1, text_label[i])
    max_value = THSvar.GetMaximum('NoStack')
    hist[list(samples)[0]].GetYaxis().SetRangeUser(0, max_value * 1.4)

    AtlasLabel.ATLASLabel(0.185, 0.885, "Internal")
    AtlasLabel.myText(0.185, 0.885-0.063, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    AtlasLabel.myText(0.185,0.885-0.063*2, region.replace("zero_", "0-").replace("one_", "1-").replace("two_", "2-") if not category_name else category_name)

    l = Leg(hist, 0.52, 0.93)
    l.Draw("SAME")

    C1.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.4)
    pad2.SetTopMargin(0)
    pad2.SetLeftMargin(0.15)
    pad2.SetBottomMargin(0.23)
    #pad2.SetGridx()
    pad2.Draw();
    pad2.cd();

    THSvar_ratio = ROOT.THStack(var + '_ratio', '')

    for sample_name in hist_ratio:
        hist_ratio[sample_name].Draw('hist E1 SAME')


    if adjustNdivision: hist_ratio[list(hist_ratio)[0]].GetXaxis().SetNdivisions(adjustNdivision, ROOT.kFALSE)

    hist_ratio[list(hist_ratio)[0]].GetXaxis().SetTitle("%s" % (varname))
    #hist_ratio[list(hist_ratio)[0]].GetXaxis().SetTitleOffset(1.5)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetNdivisions(505)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetTitle("#frac{sample}{%s}"%(denominator) if not ratio_title else ratio_title)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetTitleOffset(1.8)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetTitleSize(26)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetTitleFont(43)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetLabelSize(26)
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetLabelFont(43)
    hist_ratio[list(hist_ratio)[0]].GetXaxis().SetTitleSize(26)
    hist_ratio[list(hist_ratio)[0]].GetXaxis().SetTitleFont(43)
    hist_ratio[list(hist_ratio)[0]].GetXaxis().SetTitleOffset(2.8)
    hist_ratio[list(hist_ratio)[0]].GetXaxis().SetLabelSize(26)
    hist_ratio[list(hist_ratio)[0]].GetXaxis().SetLabelFont(43)
    if discrete:
        for i in range(nbin):
            hist_ratio[list(hist_ratio)[0]].GetXaxis().SetBinLabel(i+1, '%d' % minvar+i)
    if text_label:
         for i in range(nbin):
            hist_ratio[list(hist_ratio)[0]].GetXaxis().SetBinLabel(i+1, text_label[i])
    hist_ratio[list(hist_ratio)[0]].GetYaxis().SetRangeUser(0.7, 1.3)

    if not category:
        if not os.path.isdir('Variables_ratio/%s' % (region)):
            os.makedirs('Variables_ratio/%s' % (region))
        C1.Print("Variables_ratio/%s/%s.pdf" % (region, var.replace('/','_')))

    else:

        if not os.path.isdir('%s_ratio/%s' % (categorization, category_name)):
            os.makedirs('%s_ratio/%s' % (categorization, category_name))
        C1.Print("%s_ratio/%s/%s.pdf" % (categorization, category_name, var.replace('/','_')))

    for sample_name in samples:
        hist[sample_name].Delete()

    return

def plotVars(region, category=None, category_name=None, categorization='Event_XGB_16_Category'):

    if not category:

        samples = OrderedDict([
            ('Data sideband', {'process': ['data_sid'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kBlack}),
            # ('Data center', {'process': ['data_cen'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kBlack, 'width': 3}),
            ('VBF H#rightarrow#mu#mu simulation', {'process': ['VBF'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kRed, 'width': 3}),
            ('ggF H#rightarrow#mu#mu simulation', {'process': ['ggF'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kOrange+2, 'width': 3}),
            ('MC bkg sideband', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130))', 'color': ROOT.kBlue}),
            ('MC bkg center', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kBlue, 'width': 3})
        ])

        region_code = {'zero_jet': 0, 'one_jet': 1, 'two_jet': 2}

        bdt_boundary = {'zero_jet': [0.21, 0.53, 0.81], 'one_jet': [0.36, 0.67, 0.88], 'two_jet': [0.16, 0.42, 0.65], 'VBF': [0.62, 0.75, 0.85, 0.93]}

        # plot_var(samples, region,'ClassOut_XGB_Higgs', 'O_{ggF}^{(%d)}' % region_code[region], 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
        # if region == 'two_jet': plot_var(samples, region,'ClassOut_XGB_VBF', 'O_{VBF}', 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'ClassOut_XGB_QG_Higgs', 'O_{ggF}^{(%d)}' % region_code[region], 20, 0, 1, category=category, category_name=category_name, categorization=categorization, region_mark=bdt_boundary[region])
        if region == 'two_jet': plot_var(samples, region,'ClassOut_XGB_QG_VBF', 'O_{VBF}', 20, 0, 1, category=category, category_name=category_name, categorization=categorization, region_mark=bdt_boundary['VBF'], arrow={'ggF': [0.55, 0.3]})

    samples = OrderedDict([
        ('Data sideband', {'process': ['data_sid'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=180)', 'color': ROOT.kBlack}),
        # ('Data center', {'process': ['data_cen'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kBlack, 'width': 3}),
        ('H#rightarrow#mu#mu simulation', {'process': ['VBF', 'VH', 'ggF', 'ttH'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kRed, 'width': 3}),
        ('MC bkg sideband', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130))', 'color': ROOT.kBlue}),
        ('MC bkg center', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kBlue, 'width': 3})
    ])

    # plot_var(samples, region, 'm_mumu', 'm_{#mu#mu} [GeV]', 14, 110, 180, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'Z_PT_OnlyNearFsr', 'p_{T}^{#mu#mu} [GeV]', 20, 0, 200 if 'two_jet' in region else 100, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'Z_Y_OnlyNearFsr', 'y_{#mu#mu}', 30, -3, 3, category=category, category_name=category_name, categorization=categorization)
    # plot_var(samples, region, 'fabs(Muons_CosThetaStar)', '|cos#theta*|', 20, 0, 1, category=category, category_name=category_name, categorization=categorization)
    plot_var(samples, region, 'Muons_CosThetaStar', 'cos#theta*', 20, -1, 1, category=category, category_name=category_name, categorization=categorization)

    if region != 'zero_jet':
        plot_var(samples, region, 'Jets_PT_Lead', 'p_{T}^{j_{#scale[1]{1}}} [GeV]', 20, 0, 300, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Eta_Lead', '#eta_{j_{1}}', 20, -5, 5, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'DeltaPhi_mumuj1', '#Delta#phi_{#mu#mu,j_{1}}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_QGscore_Lead', 'N_{tracks}^{j_{1}}', 20, 0, 20, category=category, category_name=category_name, categorization=categorization)
        # plot_var(samples, region, 'Jets_QGflag_Lead', 'QG_{j_{1}}', 3, -1, 2, text_label=['Untagged', 'Quark', 'Gluon'], category=category, category_name=category_name, categorization=categorization)

    if region == 'two_jet':
        plot_var(samples, region, 'Jets_PT_Sub', 'p_{T}^{j_{#scale[1]{2}}} [GeV]', 20, 0, 200, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Eta_Sub', '#eta_{j_{2}}', 20, -4.5, 4.5, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'DeltaPhi_mumuj2', '#Delta#phi_{#mu#mu,j_{2}}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_QGscore_Sub', 'N_{tracks}^{j_{2}}', 20, 0, 20, category=category, category_name=category_name, categorization=categorization)
        # plot_var(samples, region, 'Jets_QGflag_Sub', 'QG_{j_{2}}', 3, -1, 2, text_label=['Untagged', 'Quark', 'Gluon'], category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_PT_jj', 'p_{T}^{jj} [GeV]', 20, 0, 400, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Y_jj', 'y_{jj}', 40, -5, 5, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'DeltaPhi_mumujj', '#Delta#phi_{#mu#mu,jj}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region, 'Jets_Minv_jj', 'm_{jj} [GeV]' , 20, 0, 1000, adjustNdivision = 40205, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region,'Event_MET', 'E_{T}^{miss} [GeV]', 20, 0, 150 if 'two_jet' in region else 80, category=category, category_name=category_name, categorization=categorization)
        # plot_var(samples, region,'DeltaPhi_mumuMET', '#Delta#phi_{#mu#mu,MET}', 20, -3.14, 3.14, category=category, category_name=category_name, categorization=categorization)
        plot_var(samples, region,'Event_Ht-Muons_PT_Lead-Muons_PT_Sub', 'H_{T} [GeV]', 20, 0, 400, category=category, category_name=category_name, categorization=categorization)


def plot_m_mumu(categorization='Event_Paper_Category'):

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
            17: 'VH4Lep-High',
            18: 'VH3Lep-High',
            19: 'VH3Lep-Low',
            20: 'ttH',
        }

    samples = OrderedDict([
        ('Data sideband', {'process': ['data_sid'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=160)', 'color': ROOT.kBlack}),
        ('Higgs signal', {'process': ['VBF', 'VH', 'ggF', 'ttH'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=160)', 'color': ROOT.kRed, 'width': 3}),
        ('MC bkg', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=160)', 'color': ROOT.kBlue}),
    ])

    for category in categories:

        if category <= 8: region = 'two_jet'
        elif category <= 12: region = 'one_jet'
        elif category <= 16: region = 'zero_jet'
        else: region = 'VH_ttH'

        plot_var(samples, region, 'm_mumu', 'm_{#mu#mu} [GeV]', 10, 110, 160, category=category, category_name=categories[category], categorization=categorization)

def plotVars_ratio(region, category=None, category_name=None, categorization='Event_XGB_16_Category'):

    if not category:

        samples = OrderedDict([
            # ('Data sideband', {'process': ['data_sid'], 'selection': 'weight*(m_mumu>=110&&m_mumu<=180)', 'color': ROOT.kBlack, 'ratio': 'Data center'}),
            # ('Data center', {'process': ['data_cen'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kBlue+3, 'style': 2, 'ratio': 'numerator'}),
            ('VBF', {'process': ['VBF'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kRed, 'width': 3}),
            ('ggF', {'process': ['ggF'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kOrange+2, 'width': 3}),
            ('MC bkg sideband', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130))', 'color': ROOT.kBlue, 'ratio': 'MC bkg center'}),
            ('MC bkg center', {'process': ['Z', 'ttbar', 'diboson', 'stop'], 'selection': 'weight*(m_mumu>=120&&m_mumu<=130)', 'color': ROOT.kGreen+2, 'ratio': 'numerator'})
        ])

        region_code = {'zero_jet': 0, 'one_jet': 1, 'two_jet': 2}

        plot_var_ratio(samples, region,'ClassOut_XGB_QG_Higgs', 'O_{ggF}^{(%d)}' % region_code[region], 20, 0, 1, category=category, category_name=category_name, categorization=categorization, ratio_title="#frac{center}{sideband}")
        if region == 'two_jet': plot_var_ratio(samples, region,'ClassOut_XGB_QG_VBF', 'O_{VBF}', 20, 0, 1, category=category, category_name=category_name, categorization=categorization, ratio_title="#frac{center}{sideband}")

def main():

    args=getArgs()

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.SetDefaultSumw2(1)

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

            plotVars(region, category=category, category_name=categories[category], categorization=args.categorization)

    elif args.mode == 'm_mumu':

        plot_m_mumu(categorization=args.categorization)

    elif args.mode == 'BDT_sideband_center_ratio':

        regions = ['zero_jet', 'one_jet', 'two_jet']

        for region in regions:
        
            plotVars_ratio(region)

    return

if __name__ == '__main__':
    main()
