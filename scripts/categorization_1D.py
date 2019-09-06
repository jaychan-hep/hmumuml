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
from math import sqrt, log
import math
import sys
import time
import json

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'zero_jet', help = 'Region to process')
    parser.add_argument('-n', '--nscan', type = int, default = 100, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 5, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('--floatB', action = 'store_true', default = False, help = 'Floting last boundary')
    parser.add_argument('-f', '--nfold', type = int, default = 1, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = 10, help='early stopping rounds.')
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
    start_time = time.time()

    j1max, j2max, j3max, j4max, j5max, j6max, j7max, j8max, j9max, j10max, j11max, j12max, j13max, j14max, j15max, smax = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    j8stop = 0
    for j8 in range(1, nscan+1 if nbin >= 9 else 2):

        elapsed_time = time.time() - start_time
        tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print '%s |%s%s| %d/100'%(tt, '>'*(j8*100/nscan), ' '*(100-j8*100/nscan), j8*100/nscan) # progress bar
        sys.stdout.write("\033[F") # Cursor up one line

        j1max_, j2max_, j3max_, j4max_, j5max_, j6max_, j7max_, smax_1_8 = 0, 0, 0, 0, 0, 0, 0, 0
        j4stop = 0
        for j4 in range(1, j8+1 if nbin >= 13 else 2):

            j1max__, j2max__, j3max__, j5max__, j6max__, j7max__, smax_1_4, smax_5_8 = 0, 0, 0, 0, 0, 0, 0, 0
            j2stop = 0
            for j2 in range(1, j4+1 if nbin >= 15 else 2):

                j1max___, j3max___, smax_1_2, smax_3_4 = 0, 0, 0, 0
                j1stop = 0
                for j1 in range(1, j2+1 if nbin >= 16 else 2):
                    if not (floatB and nbin == 16):
                        nsig1, dsig1 = hist_integral(h_sig, 1, j1-1) 
                        nbkg1, dbkg1 = hist_integral(h_bkg, 1, j1-1) 
                        if nbin >= 16 and nbkg1 < minN: continue
                        s1, u1 = calc_sig(nsig1, nbkg1, dsig1, dbkg1)
                    else:
                        s1, u1 = 0, 0
                    if not (floatB and nbin == 15):
                        nsig2, dsig2 = hist_integral(h_sig, j1, j2-1)
                        nbkg2, dbkg2 = hist_integral(h_bkg, j1, j2-1)
                        if nbin >= 15 and nbkg2 < minN: break
                        s2, u2 = calc_sig(nsig2, nbkg2, dsig2, dbkg2)
                    else:
                        s2, u2 = 0, 0
                    s_1_2, _ = sum_z([[s1,u1],[s2,u2]])
                    if s_1_2 >= smax_1_2:
                        j1stop = 0
                        smax_1_2 = s_1_2
                        j1max___ = j1
                    else:
                        j1stop += 1
                    if j1stop == earlystop: break

                if j1max___ == 0: continue

                j3stop = 0
                for j3 in range(j2, j4+1 if nbin >= 14 else 2):
                    if not (floatB and nbin == 14):
                        nsig3, dsig3 = hist_integral(h_sig, j2, j3-1)
                        nbkg3, dbkg3 = hist_integral(h_bkg, j2, j3-1)
                        if nbin >= 14 and nbkg3 < minN: continue
                        s3, u3 = calc_sig(nsig3, nbkg3, dsig3, dbkg3)
                    else:
                        s3, u3 = 0, 0
                    if not (floatB and nbin == 13):
                        nsig4, dsig4 = hist_integral(h_sig, j3, j4-1)
                        nbkg4, dbkg4 = hist_integral(h_bkg, j3, j4-1)
                        if nbin >= 13 and nbkg4 < minN: break
                        s4, u4 = calc_sig(nsig4, nbkg4, dsig4, dbkg4)
                    else:
                        s4, u4 = 0, 0
                    s_3_4, _ = sum_z([[s3,u3],[s4,u4]])
                    if s_3_4 >= smax_3_4:
                        j3stop = 0
                        smax_3_4 = s_3_4
                        j3max___ = j3
                    else:
                        j3stop += 1
                    if j3stop == earlystop: break

                if j3max___ == 0: break

                s_1_4 = sqrt(smax_1_2**2 + smax_3_4**2)
                if s_1_4 >= smax_1_4:
                    j2stop = 0
                    smax_1_4 = s_1_4
                    j1max__, j2max__, j3max__ = j1max___, j2, j3max___
                else:
                    j2stop += 1
                if j2stop == earlystop: break

            if j2max__ == 0: continue

            j6stop = 0
            for j6 in range(j4, j8+1 if nbin >= 11 else 2):

                j5max___, j7max___, smax_5_6, smax_7_8 = 0, 0, 0, 0
                j5stop = 0
                for j5 in range(j4, j6+1 if nbin >= 12 else 2):
                    if not (floatB and nbin == 12):
                        nsig5, dsig5 = hist_integral(h_sig, j4, j5-1)
                        nbkg5, dbkg5 = hist_integral(h_bkg, j4, j5-1)
                        if nbin >= 12 and nbkg5 < minN: continue
                        s5, u5 = calc_sig(nsig5, nbkg5, dsig5, dbkg5)
                    else:
                        s5, u5 = 0, 0
                    if not (floatB and nbin == 11):
                        nsig6, dsig6 = hist_integral(h_sig, j5, j6-1)
                        nbkg6, dbkg6 = hist_integral(h_bkg, j5, j6-1)
                        if nbin >= 11 and nbkg6 < minN: break
                        s6, u6 = calc_sig(nsig6, nbkg6, dsig6, dbkg6)
                    else:
                        s6, u6 = 0, 0
                    s_5_6, _ = sum_z([[s5,u5],[s6,u6]])
                    if s_5_6 >= smax_5_6:
                        j5stop = 0
                        smax_5_6 = s_5_6
                        j5max___ = j5
                    else:
                        j5stop += 1
                    if j5stop == earlystop: break

                if j5max___ == 0: continue

                j7stop = 0
                for j7 in range(j6, j8+1 if nbin >= 10 else 2):
                    if not (floatB and nbin == 10):
                        nsig7, dsig7 = hist_integral(h_sig, j6, j7-1)
                        nbkg7, dbkg7 = hist_integral(h_bkg, j6, j7-1)
                        if nbin >= 10 and nbkg7 < minN: continue
                        s7, u7 = calc_sig(nsig7, nbkg7, dsig7, dbkg7)
                    else:
                        s7, u7 = 0, 0
                    if not (floatB and nbin == 9):
                        nsig8, dsig8 = hist_integral(h_sig, j7, j8-1)
                        nbkg8, dbkg8 = hist_integral(h_bkg, j7, j8-1)
                        if nbin >= 9 and nbkg8 < minN: break
                        s8, u8 = calc_sig(nsig8, nbkg8, dsig8, dbkg8)
                    else:
                        s8, u8 = 0, 0
                    s_7_8, _ = sum_z([[s7,u7],[s8,u8]])
                    if s_7_8 >= smax_7_8:
                        j7stop = 0
                        smax_7_8 = s_7_8
                        j7max___ = j7
                    else:
                        j7stop += 1
                    if j7stop == earlystop: break

                if j7max___ == 0: break

                s_5_8 = sqrt(smax_5_6**2 + smax_7_8**2)
                if s_5_8 >= smax_5_8:
                    j6stop = 0
                    smax_5_8 = s_5_8
                    j5max__, j6max__, j7max__ = j5max___, j6, j7max___
                else:
                    j6stop += 1
                if j6stop == earlystop: break

            if j6max__ == 0: break

            s_1_8 = sqrt(smax_1_4**2 + smax_5_8**2)
            if s_1_8 >= smax_1_8:
                j4stop = 0
                smax_1_8 = s_1_8
                j1max_, j2max_, j3max_, j4max_, j5max_, j6max_, j7max_ = j1max__, j2max__, j3max__, j4, j5max__, j6max__, j7max__
            else:
                j4stop += 1
            if j4stop == earlystop: break

        if j4max_ == 0: continue
        ################################

        j9max_, j10max_, j11max_, j12max_, j13max_, j14max_, j15max_, smax_9_16 = 0, 0, 0, 0, 0, 0, 0, 0
        j12stop = 0
        for j12 in range(j8, nscan+1 if nbin >= 5 else 2):

            j9max__, j10max__, j11max__, j13max__, j14max__, j15max__, smax_9_12, smax_13_16 = 0, 0, 0, 0, 0, 0, 0, 0
            j10stop = 0
            for j10 in range(j8, j12+1 if nbin >= 7 else 2):

                j9max___, j11max___, smax_9_10, smax_11_12 = 0, 0, 0, 0    
                j9stop = 0
                for j9 in range(j8, j10+1 if nbin >= 8 else 2):
                    if not (floatB and nbin == 8):
                        nsig9, dsig9 = hist_integral(h_sig, j8, j9-1)
                        nbkg9, dbkg9 = hist_integral(h_bkg, j8, j9-1)
                        if nbin >= 8 and nbkg9 < minN: continue
                        s9, u9 = calc_sig(nsig9, nbkg9, dsig9, dbkg9)
                    else:
                        s9, u9 = 0, 0
                    if not (floatB and nbin == 7):
                        nsig10, dsig10 = hist_integral(h_sig, j9, j10-1)
                        nbkg10, dbkg10 = hist_integral(h_bkg, j9, j10-1)
                        if nbin >= 7 and nbkg10 < minN: break
                        s10, u10 = calc_sig(nsig10, nbkg10, dsig10, dbkg10)
                    else:
                        s10, u10 = 0, 0
                    s_9_10, _ = sum_z([[s9,u9],[s10,u10]])
                    if s_9_10 >= smax_9_10:
                        j9stop = 0
                        smax_9_10 = s_9_10
                        j9max___ = j9
                    else:
                        j9stop += 1
                    if j9stop == earlystop: break

                if j9max___ == 0: continue

                j11stop = 0
                for j11 in range(j10, j12+1 if nbin >= 6 else 2):
                    if not (floatB and nbin == 6):
                        nsig11, dsig11 = hist_integral(h_sig, j10, j11-1)
                        nbkg11, dbkg11 = hist_integral(h_bkg, j10, j11-1)
                        if nbin >= 6 and nbkg11 < minN: continue
                        s11, u11 = calc_sig(nsig11, nbkg11, dsig11, dbkg11)
                    else:
                        s11, u11 = 0, 0
                    if not (floatB and nbin == 5):
                        nsig12, dsig12 = hist_integral(h_sig, j11, j12-1)
                        nbkg12, dbkg12 = hist_integral(h_bkg, j11, j12-1)
                        if nbin >= 5 and nbkg12 < minN: break
                        s12, u12 = calc_sig(nsig12, nbkg12, dsig12, dbkg12)
                    else:
                        s12, u12 = 0, 0
                    s_11_12, _ = sum_z([[s11,u11],[s12,u12]])
                    if s_11_12 >= smax_11_12:
                        j11stop = 0
                        smax_11_12 = s_11_12
                        j11max___ = j11
                    else:
                        j11stop += 1
                    if j11stop == earlystop: break

                if j11max___ ==0: break

                s_9_12 = sqrt(smax_9_10**2 + smax_11_12**2)
                if s_9_12 >= smax_9_12:
                    j10stop = 0
                    smax_9_12 = s_9_12
                    j9max__, j10max__, j11max__ = j9max___, j10, j11max___
                else:
                    j10stop += 1
                if j10stop == earlystop: break

            if j10max__ == 0: continue

            j14stop = 0
            for j14 in range(j12, nscan+1 if nbin >= 3 else 2):

                j13max___, j15max___, smax_13_14, smax_15_16 = 0, 0, 0, 0
                j13stop = 0
                for j13 in range(j12, j14+1 if nbin >= 4 else 2):
                    if not (floatB and nbin == 4):
                        nsig13, dsig13 = hist_integral(h_sig, j12, j13-1)
                        nbkg13, dbkg13 = hist_integral(h_bkg, j12, j13-1)
                        if nbin >= 4 and nbkg13 < minN: continue
                        s13, u13 = calc_sig(nsig13, nbkg13, dsig13, dbkg13)
                    else:
                        s13, u13 = 0, 0
                    if not (floatB and nbin == 3):
                        nsig14, dsig14 = hist_integral(h_sig, j13, j14-1)
                        nbkg14, dbkg14 = hist_integral(h_bkg, j13, j14-1)
                        if nbin >= 3 and nbkg14 < minN: break
                        s14, u14 = calc_sig(nsig14, nbkg14, dsig14, dbkg14)
                    else:
                        s14, u14 = 0, 0
                    s_13_14, _ = sum_z([[s13,u13],[s14,u14]])
                    if s_13_14 >= smax_13_14:
                        j13stop = 0
                        smax_13_14 = s_13_14
                        j13max___ = j13
                    else:
                        j13stop += 1
                    if j13stop == earlystop: break

                if j13max___ == 0: continue

                j15stop = 0
                for j15 in range(j14, nscan+1 if nbin >= 2 else 2):
                    if not (floatB and nbin == 2):
                        nsig15, dsig15 = hist_integral(h_sig, j14, j15-1)
                        nbkg15, dbkg15 = hist_integral(h_bkg, j14, j15-1)
                        if nbin >= 2 and nbkg15 < minN: continue
                        s15, u15 = calc_sig(nsig15, nbkg15, dsig15, dbkg15)
                    else:
                        s15, u15 = 0, 0
                    nsig16, dsig16 = hist_integral(h_sig, j15, nscan)
                    nbkg16, dbkg16 = hist_integral(h_bkg, j15, nscan)
                    if nbkg16 < minN: break
                    s16, u16 = calc_sig(nsig16, nbkg16, dsig16, dbkg16)
                    s_15_16, _ = sum_z([[s15,u15],[s16,u16]])
                    if s_15_16 >= smax_15_16:
                        j15stop = 0
                        smax_15_16 = s_15_16
                        j15max___ = j15
                    else:
                        j15stop += 1
                    if j15stop == earlystop: break

                if j15max___ == 0: break

                s_13_16 = sqrt(smax_13_14**2 + smax_15_16**2)
                if s_13_16 >= smax_13_16:
                    j14stop = 0
                    smax_13_16 = s_13_16
                    j13max__, j14max__, j15max__ = j13max___, j14, j15max___
                else:
                    j14stop+=1
                if j14stop == earlystop: break

            if j14max__ == 0: break

            s_9_16 = sqrt(smax_9_12**2+smax_13_16**2)
            if s_9_16 >= smax_9_16:
                j12stop = 0
                smax_9_16 = s_9_16
                j9max_, j10max_, j11max_, j12max_, j13max_, j14max_, j15max_ = j9max__, j10max__, j11max__, j12, j13max__, j14max__, j15max__
            else:
                j12stop += 1
            if j12stop == earlystop: break

        if j12max_ == 0: break

        s = sqrt(smax_1_8**2 + smax_9_16**2)
        if s >= smax:
            j8stop = 0
            smax=s
            j1max, j2max, j3max, j4max, j5max, j6max, j7max, j8max, j9max, j10max, j11max, j12max, j13max, j14max, j15max = j1max_, j2max_, j3max_, j4max_, j5max_, j6max_, j7max_, j8, j9max_, j10max_, j11max_, j12max_, j13max_, j14max_, j15max_
        else:
            j8stop += 1
        if j8stop == earlystop: break


    elapsed_time = time.time() - start_time
    tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print '%s |%s| 100/100'%(tt, '>'*100)

    boundaries = [j1max, j2max, j3max, j4max, j5max, j6max, j7max, j8max, j9max, j10max, j11max, j12max, j13max, j14max, j15max]

    # remove overlapping boundaries
    i=0
#        while boundaries[0] == 1:
#            boundaries.pop(0)
    while i < len(boundaries) - 1:
        if boundaries[i+1] == boundaries[i]:
            boundaries.pop(i+1)
        else: i += 1
    if floatB and len(boundaries) == nbin: boundaries.pop(0)

    boundaries_values = [(i-1.)/nscan for i in boundaries]
    print '========================================================================='
    print 'Fold number %d' %fold
    print 'The maximal significance:  %f' %(smax)
    print 'Boundaries: ', boundaries_values
    print '========================================================================='

    return boundaries, boundaries_values, smax
    


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
