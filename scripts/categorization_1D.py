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
import ROOT
from argparse import ArgumentParser
from math import sqrt, log, ceil
import math
import sys
import time
import json
from tqdm import tqdm

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'zero_jet', help = 'Region to process')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-f', '--nfold', type = int, default = 1, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = -1, help='early stopping rounds.')
    return  parser.parse_args()

def calc_sig(sig, bkg,s_err,b_err):

  ntot = sig + bkg

  if(sig <= 0): return 0, 0
  if(bkg <= 0): return 0, 0

  #Counting experiment
  significance = sqrt(2*((ntot*log(ntot/bkg)) - sig))
  #significance = sig / sqrt(bkg)

  #error on significance
  numer = sqrt((log(ntot/bkg)*s_err)**2 + ((log(1+(sig/bkg)) - (sig/bkg))*b_err)**2)
  uncert = (numer/significance)

  return significance, uncert

def hist_integral(hist,i,j):
    err = ROOT.Double()
    if i>j:
        n=0
        err=0
    else: n = hist.IntegralAndError(i,j,err)
    return n, err

def sum_hist_integral(integrals,xs=1):
    err, n = 0, 0
    for integral in integrals:
        n += integral[0]
        err += integral[1]**2
    return n*xs, sqrt(err)*xs

def sum_z(zs):
    sumu=0
    sumz=0
    for i in range(len(zs)):
        sumz+=zs[i][0]**2
        sumu+=(zs[i][0]*zs[i][1])**2
    sumz=sqrt(sumz)
    sumu=sqrt(sumu)/sumz if sumz != 0 else 0
    return sumz,sumu

def gettingsig(region,sigs,bkgs, boundaries, nscan, nbin, n_fold, transform):

    f_bkg = ROOT.TFile('outputs/%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_bkg = []
    for fold in range(n_fold):
        h_bkg.append(ROOT.TH1F('h_bkg_%d'%(fold),'h_bkg_%d'%(fold),nscan,0.,1.))
        h_bkg[fold].Sumw2()
        t_bkg.Draw("bdt_score%s>>h_bkg_%d"%('_t' if transform else '', fold),"weight*1*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&eventNumber%%%d==%d)"%(n_fold,fold))

    nbkg = []
    for i in range(len(boundaries[0])):
        nbkg.append(sum_hist_integral([hist_integral(h_bkg[fold], boundaries[fold][i], boundaries[fold][i+1]-1 if i != len(boundaries[0])-1 else nscan+1) for fold in range(n_fold)]))    

    f_sig = ROOT.TFile('outputs/%s/sig.root' % (region))
    t_sig = f_sig.Get('test')

    h_sig = []
    for fold in range(n_fold):
        h_sig.append(ROOT.TH1F('h_sig_%d'%(fold),'h_sig_%d'%(fold),nscan,0.,1.))
        h_sig[fold].Sumw2()
        t_sig.Draw("bdt_score%s>>h_sig_%d"%('_t' if transform else '', fold),"weight*1*(m_mumu>=120&&m_mumu<=130&&eventNumber%%%d==%d)"%(n_fold,fold))

    nsig = []
    for i in range(len(boundaries[0])):
        nsig.append(sum_hist_integral([hist_integral(h_sig[fold], boundaries[fold][i], boundaries[fold][i+1]-1 if i != len(boundaries[0])-1 else nscan+1) for fold in range(n_fold)])
)
       
    zs = [ calc_sig(nsig[i][0], nbkg[i][0], nsig[i][1], nbkg[i][1]) for i in range(len(nsig))]
    s, u = sum_z(zs)

    print 'Signal yield: ', nsig
    print 'BKG yield: ', nbkg
    print 'Significance:  %f +- %f'%(s,abs(u))

    return s, abs(u)


class categorizer(object):

    def __init__(self, h_sig, h_bkg):

        self.h_sig = h_sig
        self.h_bkg = h_bkg

    def fit(self, bl, br, nbin, minN=5, floatB=False, earlystop=-1, pbar=False):
    
        if nbin == 1:

            if floatB: return [], 0

            nsig, dsig = hist_integral(self.h_sig, bl, br)
            nbkg, dbkg = hist_integral(self.h_bkg, bl, br)
    
            if nbkg < minN: return -1, -1
    
            z, u = calc_sig(nsig, nbkg, dsig, dbkg)
    
            return [], z
    
        elif nbin > 1:
    
            L = int(ceil(log(nbin, 2)))
            N2 = 2**(L - 1)
            N1 = nbin - N2
    
            bmax, zmax, stop = -1, -1, 0
    
            for b in (range(bl, br+1) if not pbar else tqdm(range(bl, br+1))):
    
                b1, z1 = self.fit(bl, b-1, N1, minN=minN, floatB=floatB, earlystop=earlystop)
                if b1 == -1: continue
                b2, z2 = self.fit(b, br, N2, minN=minN, earlystop=earlystop)
                if b2 == -1: break
    
                z = sqrt(z1**2 + z2**2)
                if z > zmax:
                    stop = 0
                    zmax = z
                    bmax = sorted(list(set(b1 + [b] + b2)))
                else:
                    stop += 1
                if stop == earlystop: break
    
            return bmax, zmax


def categorizing(region,sigs,bkgs,nscan, minN, transform, nbin, floatB, n_fold, fold, earlystop):

    f_sig = ROOT.TFile('outputs/%s/sig.root' % (region))
    t_sig = f_sig.Get('test')
 
    f_bkg = ROOT.TFile('outputs/%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig=ROOT.TH1F('h_sig','h_sig',nscan,0,1)
    h_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0,1)

    h_sig.Sumw2()
    h_bkg.Sumw2()

    t_sig.Draw("bdt_score%s>>h_sig"%('_t' if transform else ''), "weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_bkg.Draw("bdt_score%s>>h_bkg"%('_t' if transform else ''), "weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))

    print 'INFO: scanning all of the possibilities...'
    cgz = categorizer(h_sig, h_bkg)
    bmax, zmax = cgz.fit(1, nscan, nbin, minN=minN, floatB=floatB, earlystop=earlystop, pbar=True)

    boundaries = bmax if floatB else ([1] + bmax)
    boundaries_values = [(i-1.)/nscan for i in boundaries]
    print '========================================================================='
    print 'Fold number %d' %fold
    print 'The maximal significance:  %f' %(zmax)
    print 'Boundaries: ', boundaries_values
    print '========================================================================='

    return boundaries, boundaries_values, zmax
    


def main():

    ROOT.gROOT.SetBatch(True)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data_sid']

    if args.floatB and args.nbin == 16:
        print 'ERROR: With floatB option, the maximun nbin is 15!!'
        quit()


    region = args.region

    nscan=args.nscan

    CrossSections={}

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
    print 'Averaged significance: ', smax

    s, u = gettingsig(region,sigs,bkgs, boundaries, nscan, args.nbin, n_fold, args.transform)

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
        print 'INFO: Creating output folder: "significances/%s"'%region
        os.makedirs("significances/%s"%region)

    with open('significances/%s/%d.json' % (region, args.nbin), 'w') as json_file:
        json.dump(outs, json_file)

    return



if __name__ == '__main__':
    main()
