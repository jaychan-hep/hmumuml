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
import os
from ROOT import TFile, TH1F, TH1, gROOT
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from root_pandas import *
import json
from categorizer import *
from pdb import set_trace

pd.options.mode.chained_assignment = None

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'zero_jet', help = 'Region to process')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 4, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-f', '--nfold', type = int, default = 1, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = -1, help='early stopping rounds.')
    return  parser.parse_args()

def gettingsig(region, boundaries, transform):

    yields = pd.DataFrame({'sig': [0.]*len(boundaries[0]),
                          'sig_err': [0.]*len(boundaries[0]),
                          'bkg': [0.]*len(boundaries[0]),
                          'bkg_err': [0.]*len(boundaries[0]),
                          'VBF': [0.]*len(boundaries[0]),
                          'VBF_err': [0.]*len(boundaries[0]),
                          'ggF': [0.]*len(boundaries[0]),
                          'ggF_err': [0.]*len(boundaries[0])})

    for category in ['sig', 'bkg', 'VBF', 'ggF']:

        for data in tqdm(read_root('outputs/%s/%s.root' % (region, category), key='test', columns=['bdt_score%s' % ('_t' if transform else ''), 'm_mumu', 'weight', 'eventNumber'], chunksize=500000), desc='Loading %s' % category, ncols=70):
    
            if category in ['sig', 'VBF', 'ggF']:
                data = data[(data.m_mumu >= 120) & (data.m_mumu <= 130)]
                data['w'] = data.weight
    
            elif category == 'bkg':
                data = data[(data.m_mumu >= 110) & (data.m_mumu <= 180) & ((data.m_mumu < 120) | (data.m_mumu > 130))]
                data['w'] = data.weight * 0.2723
    
            for i in range(len(boundaries)):
    
                data_s = data[data.eventNumber % len(boundaries) == i]
    
                for j in range(len(boundaries[i])):
    
                    data_ss = data_s[data_s['bdt_score%s' % ('_t' if transform else '')] >= boundaries[i][j]]
                    if j != len(boundaries[i]) - 1: data_ss = data_ss[data_ss['bdt_score%s' % ('_t' if transform else '')] < boundaries[i][j+1]]
    
                    yields[category][j] += data_ss.w.sum()
                    yields[category+'_err'][j] = np.sqrt(yields[category+'_err'][j]**2 + (data_ss.w**2).sum())

    zs = calc_sig(yields.sig, yields.bkg, yields.sig_err, yields.bkg_err)
    yields['z'] = zs[0]
    yields['u'] = zs[1]
    yields['VBF purity [%]'] = yields['VBF']/yields['sig']*100
    yields['s/b'] = yields['sig']/yields['bkg']

    print(yields)

    z = np.sqrt((yields['z']**2).sum())
    u = np.sqrt((yields['z']**2 * yields['u']**2).sum())/z

    print(f'Significance:  {z:.4f} +/- {abs(u):.4f}')

    return z, abs(u)

def categorizing(region,sigs,bkgs,nscan, minN, transform, nbin, floatB, n_fold, fold, earlystop):

    f_sig = TFile('outputs/%s/sig.root' % (region))
    t_sig = f_sig.Get('test')
 
    f_bkg = TFile('outputs/%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig = TH1F('h_sig','h_sig',nscan,0,1)
    h_bkg = TH1F('h_bkg','h_bkg',nscan,0,1)

    t_sig.Draw("bdt_score%s>>h_sig"%('_t' if transform else ''), "weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_bkg.Draw("bdt_score%s>>h_bkg"%('_t' if transform else ''), "weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    print('INFO: scanning all of the possibilities...')
    cgz = categorizer(h_sig, h_bkg)
    #cgz.smooth(60, nscan)  #uncomment this line to fit a function to the BDT distribution. Usage: categorizer.smooth(left_bin_to_fit, right_bin_to_fit, SorB='S' (for signal) or 'B' (for bkg), function='Epoly2', printMessage=False (switch to "True" to print message))
    bmax, zmax = cgz.fit(1, nscan, nbin, minN=minN, floatB=floatB, earlystop=earlystop, pbar=True)

    boundaries = bmax
    boundaries_values = [(i-1.)/nscan for i in boundaries]
    print('=========================================================================')
    print(f'Fold number {fold}')
    print(f'The maximal significance:  {zmax}')
    print('Boundaries: ', boundaries_values)
    print('=========================================================================')

    return boundaries, boundaries_values, zmax
    


def main():

    gROOT.SetBatch(True)
    TH1.SetDefaultSumw2(1)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data_sid']

    if args.floatB and args.nbin == 16:
        print('ERROR: With floatB option, the maximun nbin is 15!!')
        quit()


    region = args.region

    nscan=args.nscan

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


    n_fold = args.nfold
    boundaries=[]
    boundaries_values=[]
    smaxs = []
    for j in range(n_fold):
        bound, bound_value, smax = categorizing(region, sigs, bkgs, nscan, args.minN, args.transform, args.nbin if not args.floatB else args.nbin + 1, args.floatB, n_fold, j, args.earlystop)
        boundaries.append(bound)
        boundaries_values.append(bound_value)
        smaxs.append(smax)

    smax = sum(smaxs)/n_fold
    print('Averaged significance: ', smax)

    s, u = gettingsig(region, boundaries_values, args.transform)

    outs={}
    outs['boundaries'] = boundaries
    outs['boundaries_values'] = boundaries_values
    outs['smax'] = smax
    outs['significance'] = s
    outs['Delta_significance'] = u
    outs['nscan'] = nscan
    outs['transform'] = args.transform
    outs['floatB'] = args.floatB
    outs['nfold'] = n_fold
    outs['minN'] = args.minN
    outs['fine_tuned'] = False

    if not os.path.isdir('significances/%s'%region):
        print(f'INFO: Creating output folder: "significances/{region}"')
        os.makedirs("significances/%s"%region)

    with open('significances/%s/%d.json' % (region, args.nbin), 'w') as json_file:
        json.dump(outs, json_file)

    return



if __name__ == '__main__':
    main()
