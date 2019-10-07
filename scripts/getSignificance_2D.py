#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  4.22.2019
#
#
#
#
#
import os
import ROOT
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from root_pandas import *
from tqdm import tqdm
import json
from categorizer import *

pd.options.mode.chained_assignment = None

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'two_jet', help = 'Region to process')
    parser.add_argument('-b', '--nbin', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    return  parser.parse_args()

def gettingsig(region, boundaries_VBF, boundaries, transform):

    nbin = len(boundaries[0]) + len(boundaries_VBF[0])

    yields = pd.DataFrame({'sig': [0.]*nbin,
                          'sig_err': [0.]*nbin,
                          'bkg': [0.]*nbin,
                          'bkg_err': [0.]*nbin})

    for category in ['sig', 'bkg']:

        for data in tqdm(read_root('outputs/%s/%s.root' % (region, category), key='test', columns=['bdt_score%s' % ('_t' if transform else ''), 'bdt_score_VBF%s' % ('_t' if transform else ''), 'm_mumu', 'weight', 'eventNumber'], chunksize=500000), desc='Loading signal'):

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

def main():

    ROOT.gROOT.SetBatch(True)

    args=getArgs()

    region = args.region

    define_arguement = False

    print 'Number of points: %d%%'%(len(os.listdir('cat_opt/%s/%d_%d'%(region, args.vbf, args.nbin)))/(50.)*100)

    for txt in os.listdir('cat_opt/%s/%d_%d/'%(region, args.vbf, args.nbin)):
        #print txt
        try:
            with open('cat_opt/%s/%d_%d/%s' %(region, args.vbf, args.nbin, txt)) as f:
                data = json.load(f)
        except:
            print "WARNING: corrupted file 'cat_opt/%s/%d_%d/%s'!! Will remove it."%(region, args.vbf, args.nbin, txt)
            os.remove('cat_opt/%s/%d_%d/%s' %(region, args.vbf, args.nbin, txt))
            continue
        if not define_arguement:
            transform = data['transform']
            nfold = data['nfold']
            nscanvbf = data['nscanvbf']
            nscan = data['nscan']
            floatB = data['floatB']
            minN = data['minN']
            define_arguement = True
            smax = [0 for i in range(nfold)]
            smax_raw = [0 for i in range(nfold)]
            boundaries_VBF = [0 for i in range(nfold)]
            boundaries = [0 for i in range(nfold)]
        else:
            if nscanvbf != data['nscanvbf'] or nscan != data['nscan'] or transform != data['transform'] or nfold != data['nfold'] or floatB != data['floatB']:
                print 'ERROR: files have different arguement!! Please check the files and remove the obsolete files.'
                quit()
        for fold in range(nfold):
            if data[str(fold)]['smax'] > smax[fold]:
                smax[fold] = data[str(fold)]['smax']
                smax_raw[fold] = data[str(fold)]['smax_raw']
                boundaries_VBF[fold] = data[str(fold)]['VBF']
                boundaries[fold] = data[str(fold)]['ggF']

    print 'Use transformed scores' if transform else 'Use non-transformed scores'
    print 'Number of folds: ', nfold
    print 'Number of scans of VBF boundaries: ', nscanvbf
    print 'Number of scans of ggF boundaries: ', nscan
    print 'Averaged significance of categorization sets: ', sum(smax)/nfold
    print 'Averaged non-fit significance of categorization sets: ', sum(smax_raw)/nfold

    boundaries_VBF_values = []
    boundaries_values = []
    for fold in range(nfold):
        boundaries_VBF_value = [(i-1.)/nscanvbf for i in boundaries_VBF[fold]]
        boundaries_value = [(i-1.)/nscan for i in boundaries[fold]]
        boundaries_VBF_values.append(boundaries_VBF_value)
        boundaries_values.append(boundaries_value)
        print 'Boundaries %d: ' %fold, boundaries_VBF_value, boundaries_value
    
    s, u = gettingsig(region, boundaries_VBF_values, boundaries_values, transform)

    outs={}
    outs['boundaries'] = boundaries
    outs['boundaries_values'] = boundaries_values
    outs['boundaries_VBF'] = boundaries_VBF
    outs['boundaries_VBF_values'] = boundaries_VBF_values
    outs['smax'] = sum(smax)/nfold
    outs['smax_raw'] = sum(smax_raw)/nfold
    outs['significance'] = s
    outs['Delta_significance'] = u
    outs['nscan'] = nscan
    outs['nscanvbf'] = nscanvbf
    outs['transform'] = transform
    outs['floatB'] = floatB
    outs['nfold'] = nfold
    outs['minN'] = minN
    outs['fine_tuned'] = False

    if not os.path.isdir('significances/%s'%region):
        print 'INFO: Creating output folder: "significances/%s"'%region
        os.makedirs("significances/%s"%region)

    with open('significances/%s/%d_%d.json' % (region, args.vbf, args.nbin), 'w') as json_file:
        json.dump(outs, json_file)



    return



if __name__ == '__main__':
    main()
