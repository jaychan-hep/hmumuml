#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.22.2018
#
#  Modified from ttHyy XGBosst anaylsis codes
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

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'two_jet', help = 'Region to process')
    parser.add_argument('-n', '--nscan', type = int, default = 1000, help='number of scan.')
    parser.add_argument('-b', '--nbin', type = int, default = 10, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('--minN', type = float, default = 5, help = 'minimum number of events in mass window')
    #parser.add_argument('--val', action = 'store_true', default = False, help = 'use validation samples for categroization')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument('-f', '--nfold', type = int, default = 4, help='number of folds.')
    parser.add_argument('-e', '--earlystop', type = int, default = 10, help='early stopping rounds.')
    parser.add_argument('-v', '--var', action = 'store', choices = ['Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', 'Muons_Minv_MuMu_Sigma'], default = 'Muons_Minv_MuMu_Sigma', help = 'reso var')
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

def sum_z(zs,us):
    sumu=0
    sumz=0
    for i in range(len(zs)):
        sumz+=zs[i]**2
        if us != 0: sumu+=(zs[i]*us[i])**2
    sumz=sqrt(sumz)
    if us != 0: sumu=sqrt(sumu)/sumz if sumz != 0 else 0
    return sumz,sumu

def gettingsig(region, sigs, bkgs, b_resoss, nscan, nbin, n_fold, var):

    print '=================================================='
    print '            EVALUATION IN TEST SET'
    print '=================================================='


    #if not os.path.isdir('significances'):
    #    print 'INFO: Creating output folder: \"significances\"'
    #    os.makedirs("significances")

    #txt=open('significances/%s.txt' % (region), 'w')
    #for fold in range(n_fold):
    #    txt.write('-----------------------\n')
    #    txt.write('Fold number %d\n'%fold)
    #    txt.write('-----------------------\n')
    #    txt.write('#nscan: %d\n'%nscan)
    #    txt.write('#hmaxD: %d\n'%hmaxD[fold])
    #    txt.write('#First_Bondary: %f\n'%((hmaxD[fold]-1.)/nscan))
    #    txt.write('_______________________________________________________________________________________________________________________________________\n\n')

    no_reso = []
    for cat in range(1, nbin+1):
        for fold in range(n_fold):
            if b_resoss[fold][cat-1] == 0:
                no_reso.append(cat)
                break


    f_bkg = ROOT.TFile('outputs_2D/model_%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')



    #h_sig=ROOT.TH1F('h_sig','h_sig',nscan,0.01,0.03)
    #h_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0.01,0.03)
    #hreso_sig=ROOT.TH1F('hreso_sig','hreso_sig',nscan,0.01,0.03)
    #h_sigw=ROOT.TH1F('h_sigw','h_sigw',nscan,0.01,0.03)
    #hreso_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0,3)

    #h_sig.Sumw2()
    #h_bkg.Sumw2()
    #hreso_sig.Sumw2()
    #h_sigw.Sumw2()
    #hreso_bkg.Sumw2()

    h_bkg = []
    for fold in range(n_fold):
        h_bkg.append([])
        for cat in range(1, nbin+1):
            h_bkg[fold].append(ROOT.TH1F('h_bkg_%d_%d'%(fold, cat),'h_bkg_%d_%d'%(fold, cat),nscan,0.01,0.03))
            h_bkg[fold][cat-1].Sumw2()
            t_bkg.Draw("%s/m_mumu>>h_bkg_%d_%d"%(var, fold, cat),"weight*(235005./769468.)*((m_mumu>=110&&m_mumu<=160)&&!(m_mumu>=120&&m_mumu<=130)&&(bdt_category==%d)&&(eventNumber%%%d==%d))"%(cat, n_fold,fold))

    nbkgs, dbkgs = [], []
    nbkg1s, dbkg1s = [], []
    nbkg2s, dbkg2s = [], []
    for cat in range(1, nbin+1):
        nbkgs.append(sum_hist_integral([hist_integral(h_bkg[fold][cat-1], -1, nscan+1) for fold in range(n_fold)])[0])
        if not cat in no_reso:
            nbkg1s.append(sum_hist_integral([hist_integral(h_bkg[fold][cat-1], -1, b_resoss[fold][cat-1]-1) for fold in range(n_fold)])[0])
            nbkg2s.append(sum_hist_integral([hist_integral(h_bkg[fold][cat-1], b_resoss[fold][cat-1], nscan+1) for fold in range(n_fold)])[0])
        else:
            nbkg1s.append(0)
            nbkg2s.append(0)


    f_sig = ROOT.TFile('outputs_2D/model_%s/sig.root' % (region))
    t_sig = f_sig.Get('test')

    h_sig = []
    h_sigw = []
    hreso_sig = []
    for fold in range(n_fold):
        h_sig.append([])
        h_sigw.append([])
        hreso_sig.append([])
        for cat in range(1, nbin+1):
            h_sig[fold].append(ROOT.TH1F('h_sig_%d_%d'%(fold, cat),'h_sig_%d'%(fold),nscan,0.01,0.03))
            h_sigw[fold].append(ROOT.TH1F('h_sigw_%d_%d'%(fold, cat),'h_sigw_%d'%(fold),nscan,0.01,0.03))
            hreso_sig[fold].append(ROOT.TH1F('hreso_sig_%d_%d'%(fold, cat),'hreso_sig_%d'%(fold),nscan,0.01,0.03))
            h_sig[fold][cat-1].Sumw2()
            h_sigw[fold][cat-1].Sumw2()
            hreso_sig[fold][cat-1].Sumw2()

            t_sig.Draw("%s/m_mumu>>h_sig_%d_%d"%(var, fold, cat),"weight*((m_mumu>=120&&m_mumu<=130)&&(bdt_category==%d)&&(eventNumber%%%d==%d))"%(cat, n_fold,fold))
            t_sig.Draw("%s/m_mumu>>h_sigw_%d_%d"%(var, fold, cat),"weight*((m_mumu>=110&&m_mumu<=140)&&(bdt_category==%d)&&(eventNumber%%%d==%d))"%(cat, n_fold,fold))
            t_sig.Draw("%s/m_mumu>>hreso_sig_%d_%d"%(var, fold, cat),"%s*weight*((m_mumu>=110&&m_mumu<=140)&&(bdt_category==%d)&&(eventNumber%%%d==%d))"%(var, cat, n_fold,fold))
       
    nsigs, dsigs = [], []
    nsigsw, dsigsw = [], []
    nsig1s, dsig1s = [], []
    nsig2s, dsig2s = [], []
    nsig1sw, dsig1sw = [], []
    nsig2sw, dsig2sw = [], []
    ave_resos, r1s, r2s = [], [], []
    z0s, z1s, z2s, zs = [], [], [], []
    for cat in range(1, nbin+1):
        nsigs.append(sum_hist_integral([hist_integral(h_sig[fold][cat-1], -1, nscan+1) for fold in range(n_fold)])[0])
        nsigsw.append(sum_hist_integral([hist_integral(h_sigw[fold][cat-1], -1, nscan+1) for fold in range(n_fold)])[0])
        z0s.append(calc_sig(nsigs[cat-1], nbkgs[cat-1], 0, 0)[0])
        ave_resos.append(sum_hist_integral([hist_integral(hreso_sig[fold][cat-1], -1, nscan+1) for fold in range(n_fold)])[0]/nsigsw[cat-1])
        if not cat in no_reso:
            nsig1s.append(sum_hist_integral([hist_integral(h_sig[fold][cat-1], -1, b_resoss[fold][cat-1]-1) for fold in range(n_fold)])[0])
            nsig1sw.append(sum_hist_integral([hist_integral(h_sigw[fold][cat-1], -1, b_resoss[fold][cat-1]-1) for fold in range(n_fold)])[0])
            r1s.append((sum_hist_integral([hist_integral(hreso_sig[fold][cat-1], -1, b_resoss[fold][cat-1]-1) for fold in range(n_fold)])[0]/nsig1sw[cat-1])/(ave_resos[cat-1]))
            z1s.append(calc_sig(nsig1s[cat-1], nbkg1s[cat-1], 0, 0)[0]/sqrt(r1s[cat-1]))
            nsig2s.append(sum_hist_integral([hist_integral(h_sig[fold][cat-1], b_resoss[fold][cat-1], nscan+1) for fold in range(n_fold)])[0])
            nsig2sw.append(sum_hist_integral([hist_integral(h_sigw[fold][cat-1], b_resoss[fold][cat-1], nscan+1) for fold in range(n_fold)])[0])
            r2s.append((sum_hist_integral([hist_integral(hreso_sig[fold][cat-1], b_resoss[fold][cat-1], nscan+1) for fold in range(n_fold)])[0]/nsig2sw[cat-1])/(ave_resos[cat-1]))
            z2s.append(calc_sig(nsig2s[cat-1], nbkg2s[cat-1], 0, 0)[0]/sqrt(r2s[cat-1]))
            zs.append(sqrt(z1s[cat-1]**2+z2s[cat-1]**2))
        else:
            nsig1s.append(0)
            nsig1sw.append(0)
            r1s.append(0)
            z1s.append(0)
            nsig2s.append(0)
            nsig2sw.append(0)
            r2s.append(0)
            z2s.append(0)
            zs.append(z0s[cat-1])

        print 'BDT_%s_%d' %(region, cat)
        print 'Original Significance: ', z0s[cat-1]
        if not cat in no_reso:
            print 'r1: %f      r2: %f' %(r1s[cat-1], r2s[cat-1])
            print 'Significance W/ reso cat: ', zs[cat-1]
        else:
            print 'No resolution category!'
        print '-------------------------------------'

    z0 = sum_z(z0s, 0)[0]
    z = sum_z(zs, 0)[0]
       
    print 'Original Significance:    %f'%(z0)
    print 'Significance W/ reso cat: %f'%(z)


def categorizing(region, sigs, bkgs, nscan, minN, cat, n_fold, fold, earlystop, var):

    print '========================================================='
    print '   Start resolution categorization for %s_%d'%(region, cat)

    f_sig = ROOT.TFile('outputs_2D/model_%s/sig.root' % (region))
    t_sig = f_sig.Get('test')
 
    f_bkg = ROOT.TFile('outputs_2D/model_%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig=ROOT.TH1F('h_sig','h_sig',nscan,0.01,0.03)
    h_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0.01,0.03)
    hreso_sig=ROOT.TH1F('hreso_sig','hreso_sig',nscan,0.01,0.03)
    h_sigw=ROOT.TH1F('h_sigw','h_sigw',nscan,0.01,0.03)
    #hreso_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0,0.04)


    h_sig.Sumw2()
    h_bkg.Sumw2()
    hreso_sig.Sumw2()
    h_sigw.Sumw2()
    #hreso_bkg.Sumw2()

    t_sig.Draw("%s/m_mumu>>h_sig"%(var),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(bdt_category==%d)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.), cat, n_fold,fold))
    t_bkg.Draw("%s/m_mumu>>h_bkg"%(var),"weight*%f*(235005./769468.)*((m_mumu>=110&&m_mumu<=160)&&!(m_mumu>=120&&m_mumu<=130)&&(bdt_category==%d)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.), cat, n_fold,fold))
    t_sig.Draw("%s/m_mumu>>hreso_sig"%(var),"%s*weight*%f*((m_mumu>=110&&m_mumu<=140)&&(bdt_category==%d)&&(eventNumber%%%d!=%d))"%(var, n_fold/(n_fold-1.), cat, n_fold,fold))
    t_sig.Draw("%s/m_mumu>>h_sigw"%(var),"weight*%f*((m_mumu>=110&&m_mumu<=140)&&(bdt_category==%d)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.), cat, n_fold,fold))
    #t_bkg.Draw("%s>>hreso_bkg"%(var),"%s*weight*%f*(235005./769468.)*((m_mumu>=110&&m_mumu<=160)&&!(m_mumu>=120&&%s<=130)&&(eventNumber%%%d!=%d))"%(var, n_fold/(n_fold-1.),n_fold,fold))

    reso_max=0

    nsig, dsig = hist_integral(h_sig, 0, nscan+1)
    nsigw, dsigw = hist_integral(h_sigw, 0, nscan+1)
    nbkg, dbkg = hist_integral(h_bkg, 0, nscan+1)
    smax, umax = calc_sig(nsig, nbkg, dsig, dbkg)
    ave_reso = hreso_sig.Integral(0, nscan+1)/nsigw


    print 'Original significance: %f' %(smax)
    print 'nsig: ', nsig
    print 'wider nsig: ', nsigw
    print 'nbkg: ', nbkg
    print 'Average resolution ratio: %f' %(ave_reso/125.)
    print 'INFO: scanning all of the possibilities...'
    start_time = time.time()

    #Bstop = 0
    for reso in range(1, nscan+1):

        elapsed_time = time.time() - start_time
        tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print '%s |%s%s| %d/100'%(tt, '>'*(reso*100/nscan), ' '*(100-reso*100/nscan), reso*100/nscan)
        sys.stdout.write("\033[F") # Cursor up one line
        #time.sleep(0.001)

        nsig1, dsig1 = hist_integral(h_sig, 0, reso-1)
        nsig1w, dsig1w = hist_integral(h_sigw, 0, reso-1)
        nbkg1, dbkg1 = hist_integral(h_bkg, 0, reso-1)
        if nbkg1 < minN or nsig1w == 0: continue
        r1 = (hreso_sig.Integral(0, reso-1)/nsig1w)/ave_reso
        if r1==0: print 'What?!!'
        s1, u1 = calc_sig(nsig1, nbkg1, dsig1, dbkg1)
        s1 = s1/sqrt(r1)
        nsig2, dsig2 = hist_integral(h_sig, reso, nscan+1)
        nsig2w, dsig2w = hist_integral(h_sigw, reso, nscan+1)
        nbkg2, dbkg2 = hist_integral(h_bkg, reso, nscan+1)
        if nbkg2 < minN or nsig2w == 0: break
        r2 = (hreso_sig.Integral(reso, nscan+1)/nsig2w)/ave_reso
        s2, u2 = calc_sig(nsig2, nbkg2, dsig2, dbkg2)
        s2 = s2/sqrt(r2)
        s, u = sum_z([s1,s2],[u1,u2])
    
        if s >= smax:
            smax = s
            nsig1max = nsig1
            nsig2max = nsig2
            nsig1wmax = nsig1w
            nsig2wmax = nsig2w
            reso_max = reso

    elapsed_time = time.time() - start_time
    tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    print '%s |%s| 100/100'%(tt, '>'*100)
    print 'The maximal significance:  %f' %(smax)
    if reso_max != 0:
        print 'Resolution boundary:  %f' %((reso_max-1.)/nscan*0.02+0.01)
        print 'nsig1: ', nsig1max
        print 'nsig2: ', nsig2max
        print 'nsig1w: ', nsig1wmax
        print 'nsig2w: ', nsig2wmax
    else:
        print 'No resolution category!'

    print '========================================================================='

    return reso_max


def main():

    ROOT.gROOT.SetBatch(True)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data']




    region = args.region

    nscan=args.nscan

    CrossSections={}

    if not args.skip:
        siglist=''
        for sig in sigs:
            if os.path.isfile('outputs_2D/model_%s/%s.root'% (region,sig)): siglist+=' outputs_2D/model_%s/%s.root'% (region,sig)
        os.system("hadd -f outputs_2D/model_%s/sig.root"%(region)+siglist)

    #if not os.path.isfile('outputs/model_%s_%s/bkg.root' % (train_sig,region)):
    if not args.skip:
        bkglist=''
        for bkg in bkgs:
            if os.path.isfile('outputs_2D/model_%s/%s.root'% (region,bkg)): bkglist+=' outputs_2D/model_%s/%s.root'% (region,bkg)
        os.system("hadd -f outputs_2D/model_%s/bkg.root"%(region)+bkglist)


    n_fold = args.nfold
    nbin = args.nbin

    b_resoss = []
    value_resoss = []
    for fold in range(n_fold):
        b_resos = []
        value_resos = []
        for cat in range(1, nbin+1):
            b_reso = categorizing(region, sigs, bkgs, nscan, args.minN, cat, n_fold, fold, args.earlystop, args.var)
            b_resos.append(b_reso)
            value_resos.append((b_reso-1.)/nscan*0.02+0.01)
        b_resoss.append(b_resos)
        value_resoss.append(value_resos)

    print value_resoss

    gettingsig(region, sigs, bkgs, b_resoss, nscan, nbin, n_fold, args.var)

    return



if __name__ == '__main__':
    main()
