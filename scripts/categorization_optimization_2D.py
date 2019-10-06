#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.22.2018
#
#
#
#
#
#
import os
from ROOT import *
from argparse import ArgumentParser
import json
from categorizer import *

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'two_jet', help = 'Region to process')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('--nscanvbf', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 5, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 5, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--vb', type = int, default = 60, help='boundary between VBF and ggF.')
    parser.add_argument('--skip', action = 'store_true', default = True, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-f', '--nfold', type = int, default = 4, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = 30, help='early stopping rounds.')
    return  parser.parse_args()

def categorizing(region,sigs,bkgs,nscan, nscanvbf, minN, transform, nbin, vb, nvbf, floatB, n_fold, fold, earlystop):

    # get inputs
    f_sig = TFile('outputs/%s/sig.root' % (region))
    t_sig = f_sig.Get('test')
 
    f_bkg = TFile('outputs/%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig=TH1F('h_sig','h_sig',nscanvbf,0,1)
    h_sig_raw=TH2F('h_sig_raw','h_sig_raw',nscanvbf,0,1,nscan,0,1)
    h_sig_raw.Sumw2()
    h_sig.Sumw2()

    # filling signal histograms
    t_sig.Draw("bdt_score%s:bdt_score_VBF%s>>h_sig_raw"%('_t' if transform else '', '_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_sig.Draw("bdt_score_VBF%s>>h_sig"%('_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,n_fold,fold if n_fold != 1 else 1))

    h_bkg_raw = TH2F('h_bkg_raw', 'h_bkg_raw', nscanvbf, 0., 1., nscan, 0., 1.)
    h_bkg = TH1F('h_bkg', 'h_bkg', nscanvbf, 0., 1.)
    h_bkg_raw.Sumw2()
    h_bkg.Sumw2()

    # filling bkg histograms
    t_bkg.Draw("bdt_score%s:bdt_score_VBF%s>>h_bkg_raw"%('_t' if transform else '','_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_bkg.Draw("bdt_score_VBF%s>>h_bkg"%('_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    #================================================
    # categorization for VBF categories. Will scan all of the possibilities of the BDT boundaries and get the one that gives the largest significance

    cgz = categorizer(h_sig, h_bkg)
    h_sig.Delete()
    h_bkg.Delete()
    fitboundary = 0.5
    cgz.smooth(int(fitboundary * nscanvbf + 1), nscanvbf)
    boundaries_VBF, svmax = cgz.fit(vb, nscanvbf, nvbf, minN=minN, earlystop=earlystop, pbar=True)

    boundaries_VBF_values = [(i-1.)/nscanvbf for i in boundaries_VBF]

    print '========================================================================='
    print 'Fold number %d' %fold
    print 'The maximal VBF category significance:  %f' %(svmax)
    print 'VBF boundaries: ', boundaries_VBF_values

    #=================================================
    # categorization for general higgs categories

    h_sig=TH1F('h_sig','h_sig',nscan,0,1)
    h_sig.Sumw2()
    t_sig.Draw("bdt_score%s>>h_sig"%('_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,n_fold,fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf))

    h_bkg = TH1F('h_bkg', 'h_bkg', nscan, 0., 1.)
    h_bkg.Sumw2()
    t_bkg.Draw("bdt_score%s>>h_bkg"%('_t' if transform else ''), "weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf))

    cgz = categorizer(h_sig, h_bkg)
    fitboundary = 0.5
    cgz.smooth(int(fitboundary * nscan + 1), nscan)
    boundaries, smax = cgz.fit(1, nscan, nbin, minN=minN, floatB=floatB, earlystop=earlystop, pbar=True)

    boundaries_values = [(i-1.)/nscan for i in boundaries]
    print 'The maximal ggF category significance:  %f' %(smax)
    print 'ggF boundaries: ', boundaries_values
    smaxtot = sqrt(svmax**2+smax**2)
    print 'Total significane by fit: ', smaxtot


    # get the significance estimated from the unfitted histogram
    nbkg = []
    for i in range(len(boundaries_VBF)):
        nbkg.append(hist_integral(h_bkg_raw, boundaries_VBF[i], boundaries_VBF[i+1]-1 if i != len(boundaries_VBF)-1 else nscanvbf+1, 1, nscan+1))
    for i in range(len(boundaries)):
        nbkg.append(hist_integral(h_bkg_raw, 1, boundaries_VBF[0]-1, boundaries[i], boundaries[i+1]-1 if i != len(boundaries)-1 else nscan+1))


    nsig = []
    for i in range(len(boundaries_VBF)):
        nsig.append(hist_integral(h_sig_raw, boundaries_VBF[i], boundaries_VBF[i+1]-1 if i != len(boundaries_VBF)-1 else nscanvbf+1, 1, nscan+1))
    for i in range(len(boundaries)):
        nsig.append(hist_integral(h_sig_raw, 1, boundaries_VBF[0]-1, boundaries[i], boundaries[i+1]-1 if i != len(boundaries)-1 else nscan+1))

    zs = [ calc_sig(nsig[i][0], nbkg[i][0], nsig[i][1], nbkg[i][1]) for i in range(len(nsig))]
    s, u = sum_z(zs)

    print 'Significance from raw event yields in categorization set: ', s
    print '========================================================================='

    #raw_input('Press enter to continue')

    return boundaries_VBF, boundaries, smaxtot, s


def main():

    gROOT.SetBatch(True)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data']

    bkgs_VBF = ['data','ggF']




    region = args.region

    nscan=args.nscan
    nscanvbf = args.nscanvbf


    n_fold = args.nfold
    outs={}
    for j in range(n_fold):
        out = categorizing(region, sigs, bkgs, nscan, nscanvbf, args.minN, args.transform, args.nbin if not args.floatB else args.nbin + 1, args.vb, args.vbf, args.floatB, n_fold, j, args.earlystop)
        outs[j]={}
        outs[j]['VBF'] = out[0]
        outs[j]['ggF'] = out[1]
        outs[j]['smax'] = out[2]
        outs[j]['smax_raw'] = out[3]

    outs['transform'] = args.transform
    outs['nfold'] = n_fold
    outs['nscanvbf'] = nscanvbf
    outs['nscan'] = nscan
    outs['floatB'] = args.floatB
    outs['minN'] = args.minN
    #outs['region'] = region

    if not os.path.isdir('cat_opt/%s/%d_%d'%(region, args.vbf, args.nbin)):
        print 'INFO: Creating directory: cat_opt/%s/%d_%d' %(region, args.vbf, args.nbin)
        os.makedirs('cat_opt/%s/%d_%d'%(region, args.vbf, args.nbin))

    with open('cat_opt/%s/%d_%d/%d.json' % (region, args.vbf, args.nbin, args.vb), 'w') as json_file:
        json.dump(outs, json_file)

    return



if __name__ == '__main__':
    main()
