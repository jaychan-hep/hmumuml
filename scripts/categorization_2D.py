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
import numpy as np
import pandas as pd
from root_pandas import *
import json
from categorizer import *
from tqdm import tqdm

pd.options.mode.chained_assignment = None

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'two_jet', help = 'Region to process')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('--nscanvbf', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-f', '--nfold', type = int, default = 1, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = 30, help='early stopping rounds.')
    return  parser.parse_args()

def gettingsig(region, boundaries_VBF, boundaries, transform):

    nbin = len(boundaries[0]) + len(boundaries_VBF[0])

    yields = pd.DataFrame({'sig': [0.]*nbin,
                          'sig_err': [0.]*nbin,
                          'bkg': [0.]*nbin,
                          'bkg_err': [0.]*nbin})

    for category in ['sig', 'bkg']:

        for data in tqdm(read_root('outputs/%s/%s.root' % (region, category), key='test', columns=['bdt_score%s' % ('_t' if transform else ''), 'bdt_score_VBF%s' % ('_t' if transform else ''), 'm_mumu', 'weight', 'eventNumber'], chunksize=500000), desc='Loading %s' % category):
    
            if category == 'sig':
                data = data[(data.m_mumu >= 120) & (data.m_mumu <= 130)]
                data['w'] = data.weight
    
            elif category == 'bkg':
                data = data[(data.m_mumu >= 110) & (data.m_mumu <= 180) & ((data.m_mumu < 120) | (data.m_mumu > 130))]
                data['w'] = data.weight * 0.2723
    
            for i in range(len(boundaries)):
    
                data_s = data[data.eventNumber % len(boundaries) == i]
    
    
                data_ss = data_s[data_s['bdt_score_VBF%s' % ('_t' if transform else '')] < boundaries_VBF[i][0]]
    
                for j in range(len(boundaries[i])):
    
                    data_sss = data_ss[data_ss['bdt_score%s' % ('_t' if transform else '')] >= boundaries[i][j]]
                    if j != len(boundaries[i]) - 1: data_sss = data_sss[data_sss['bdt_score%s' % ('_t' if transform else '')] < boundaries[i][j+1]]
    
                    yields[category][j] += data_sss.w.sum()
                    yields[category+'_err'][j] = np.sqrt(yields[category+'_err'][j]**2 + (data_sss.w**2).sum())
    
    
                for j in range(len(boundaries_VBF[i])):
    
                    data_ss = data_s[data_s['bdt_score_VBF%s' % ('_t' if transform else '')] >= boundaries_VBF[i][j]]
                    if j != len(boundaries_VBF[i]) - 1: data_ss = data_ss[data_ss['bdt_score_VBF%s' % ('_t' if transform else '')] < boundaries_VBF[i][j+1]]
    
                    yields[category][j + len(boundaries[i])] += data_ss.w.sum()
                    yields[category+'_err'][j + len(boundaries[i])] = np.sqrt(yields[category+'_err'][j]**2 + (data_ss.w**2).sum())


    zs = calc_sig(yields.sig, yields.bkg, yields.sig_err, yields.bkg_err)
    yields['z'] = zs[0]
    yields['u'] = zs[1]

    print yields

    z = np.sqrt((yields['z']**2).sum())
    u = np.sqrt((yields['z']**2 * yields['u']**2).sum())/z

    print 'Significance:  %f +/- %f' % (z, abs(u))

    return z, abs(u)


def categorizing(region,sigs,bkgs,nscan, nscanvbf, minN, transform, nbin, nvbf, floatB, n_fold, fold, earlystop):

    # get inputs
    f_sig = TFile('outputs/%s/sig.root' % (region))
    t_sig = f_sig.Get('test')

    f_bkg = TFile('outputs/%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig_raw = TH2F('h_sig_raw', 'h_sig_raw', nscanvbf, 0., 1., nscan, 0., 1.) 
    h_sig=TH1F('h_sig','h_sig',nscanvbf,0,1)

    # filling signal histograms
    t_sig.Draw("bdt_score%s:bdt_score_VBF%s>>h_sig_raw"%('_t' if transform else '', '_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_sig.Draw("bdt_score_VBF%s>>h_sig"%('_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,n_fold,fold if n_fold != 1 else 1))

    h_bkg_raw = TH2F('h_bkg_raw', 'h_bkg_raw', nscanvbf, 0., 1., nscan, 0., 1.)
    h_bkg = TH1F('h_bkg', 'h_bkg', nscanvbf, 0., 1.)

    # filling bkg histograms
    t_bkg.Draw("bdt_score%s:bdt_score_VBF%s>>h_bkg_raw"%('_t' if transform else '','_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_bkg.Draw("bdt_score_VBF%s>>h_bkg"%('_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    fitboundary = 0.5
    fitboundary_g = 0.5

    cgz = categorizer(h_sig, h_bkg)
    cgz.smooth(int(fitboundary * nscanvbf + 1), nscanvbf)

    zmax, boundaries_VBF, boundaries = 0, -1, -1
    for vb in tqdm(range(45, 95)):

        boundaries_VBF, zv = cgz.fit(vb, nscanvbf, nvbf, minN=minN, earlystop=earlystop)
    
        h_sig_g = TH1F('h_sig_g','h_sig_g',nscan,0,1)
        t_sig.Draw("bdt_score%s>>h_sig_g"%('_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,n_fold,fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf))
    
        h_bkg_g = TH1F('h_bkg_g', 'h_bkg_g', nscan, 0., 1.)
        t_bkg.Draw("bdt_score%s>>h_bkg_g"%('_t' if transform else ''), "weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf))
    
        cgz_g = categorizer(h_sig_g, h_bkg_g)

        h_sig_g.Delete()
        h_bkg_g.Delete()

        cgz_g.smooth(int(fitboundary_g * nscan + 1), nscan)
        boundaries, zg = cgz_g.fit(1, nscan, nbin, minN=minN, floatB=floatB, earlystop=earlystop)
    
        z = sqrt(zv**2 + zg**2)
        if z >= zmax:
            zmax, boundaries_VBF_max, boundaries_max = z, boundaries_VBF, boundaries

    boundaries_VBF_values = [(i-1.)/nscanvbf for i in boundaries_VBF_max]
    boundaries_values = [(i-1.)/nscan for i in boundaries_max]

    print '========================================================================='
    print 'Fold number %d' %fold
    print 'VBF boundaries: ', boundaries_VBF_values
    print 'ggF boundaries: ', boundaries_values
    print 'Total significane by fit: ', zmax


    # get the significance estimated from the unfitted histogram
    nbkg = []
    for i in range(len(boundaries_VBF_max)):
        nbkg.append(hist_integral(h_bkg_raw, boundaries_VBF_max[i], boundaries_VBF_max[i+1]-1 if i != len(boundaries_VBF_max)-1 else nscanvbf+1, 1, nscan+1))
    for i in range(len(boundaries_max)):
        nbkg.append(hist_integral(h_bkg_raw, 1, boundaries_VBF_max[0]-1, boundaries_max[i], boundaries_max[i+1]-1 if i != len(boundaries_max)-1 else nscan+1))


    nsig = []
    for i in range(len(boundaries_VBF_max)):
        nsig.append(hist_integral(h_sig_raw, boundaries_VBF_max[i], boundaries_VBF_max[i+1]-1 if i != len(boundaries_VBF_max)-1 else nscanvbf+1, 1, nscan+1))
    for i in range(len(boundaries_max)):
        nsig.append(hist_integral(h_sig_raw, 1, boundaries_VBF_max[0]-1, boundaries_max[i], boundaries_max[i+1]-1 if i != len(boundaries_max)-1 else nscan+1))

    zs = [ calc_sig(nsig[i][0], nbkg[i][0], nsig[i][1], nbkg[i][1]) for i in range(len(nsig))]
    z, u = sum_z(zs)

    print 'Significance from raw event yields in categorization set: ', z
    print '========================================================================='

    #raw_input('Press enter to continue')

    return boundaries_VBF_max, boundaries_VBF_values, boundaries_max, boundaries_values, zmax, z

def hist_integral(hist,i,j,k=-999,l=-999):
    err = Double()
    if i>j or k>l:
        n=0
        err=0
    elif k == -999 and l == -999: n = hist.IntegralAndError(i,j,err)
    else: n = hist.IntegralAndError(i,j,k,l,err)
    return n, err

def main():

    gROOT.SetBatch(True)
    TH1.SetDefaultSumw2(1)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data_sid']

    region = args.region

    if not args.skip:
        siglist=''
        for sig in sigs:
            if os.path.isfile('outputs/%s/%s.root'% (region,sig)): siglist+=' outputs/%s/%s.root'% (region,sig)
        os.system("hadd -f outputs/%s/sig.root"%(region)+siglist)

    #if not os.path.isfile('outputs/model_%s_%s/bkg.root' % (train_sig,region)):
    if not args.skip:
        bkglist=''
        for bkg in bkgs:
            if os.path.isfile('outputs/%s/%s.root'% (region,bkg)): bkglist+=' outputs/%s/%s.root'% (region,bkg)
        os.system("hadd -f outputs/%s/bkg.root"%(region)+bkglist)

    nscan=args.nscan
    nscanvbf = args.nscanvbf

    n_fold = args.nfold
    boundaries=[]
    boundaries_values=[]
    boundaries_VBF=[]
    boundaries_VBF_values=[]
    smaxs = []
    smaxs_raw = []
    for j in range(n_fold):
        out = categorizing(region, sigs, bkgs, nscan, nscanvbf, args.minN, args.transform, args.nbin if not args.floatB else args.nbin + 1, args.vbf, args.floatB, n_fold, j, args.earlystop)
        boundaries_VBF.append(out[0])
        boundaries_VBF_values.append(out[1])
        boundaries.append(out[2])
        boundaries_values.append(out[3])
        smaxs.append(out[4])
        smaxs_raw.append(out[5])

    smax = sum(smaxs)/n_fold
    print 'Averaged significance: ', smax
    smax_raw = sum(smaxs_raw)/n_fold
    print 'Averaged raw significance: ', smax_raw

    s, u = gettingsig(region, boundaries_VBF_values, boundaries_values, args.transform)

    outs={}
    outs['boundaries'] = boundaries
    outs['boundaries_values'] = boundaries_values
    outs['boundaries_VBF'] = boundaries_VBF
    outs['boundaries_VBF_values'] = boundaries_VBF_values
    outs['smax'] = smax
    outs['smax_raw'] = smax_raw
    outs['significance'] = s
    outs['Delta_significance'] = u
    outs['nscan'] = nscan
    outs['nscanvbf'] = nscanvbf
    outs['transform'] = args.transform
    outs['floatB'] = args.floatB
    outs['nfold'] = n_fold
    outs['minN'] = args.minN
    outs['fine_tuned'] = False

    if not os.path.isdir('significances/%s'%region):
        print 'INFO: Creating output folder: "significances/%s"'%region
        os.makedirs("significances/%s"%region)

    with open('significances/%s/%d_%d.json' % (region, args.vbf, args.nbin), 'w') as json_file:
        json.dump(outs, json_file) 


if __name__ == '__main__':
    main()
