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

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-n','--nscan',type=int,default=100,help='number of scan.')
    parser.add_argument('--minN',type=float,default=2,help='minimum number of events in mass window')
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
        sumu+=(zs[i]*us[i])**2
    sumz=sqrt(sumz)
    sumu=sqrt(sumu)/sumz if sumz != 0 else 0
    return sumz,sumu

def gettingsig(region,sigs,bkgs, hmax,imax,jmax,kmax,lmax,mmax,nmax,nscan):

    var='m_mumu'
  
    if not os.path.isdir('significances'):
        print 'INFO: Creating output folder: \"significances\"'
        os.makedirs("significances")

    txt=open('significances/%s.txt' % (region), 'w')
    txt.write('#nscan: %d\n'%nscan)
    txt.write('#hmax: %d\n'%hmax)
    txt.write('#imax: %d\n'%imax)
    txt.write('#jmax: %d\n'%jmax)
    txt.write('#kmax: %d\n'%kmax)
    txt.write('#lmax: %d\n'%lmax)
    txt.write('#mmax: %d\n'%mmax)
    txt.write('#nmax: %d\n'%nmax)
    txt.write('#First_Bondary: %f\n'%((hmax-1.)/nscan))
    txt.write('#Second_Bondary: %f\n'%((imax-1.)/nscan))
    txt.write('#Third_Bondary: %f\n'%((jmax-1.)/nscan))
    txt.write('#Fourth_Bondary: %f\n'%((kmax-1.)/nscan))
    txt.write('#Fifth_Bondary: %f\n'%((lmax-1.)/nscan))
    txt.write('#Sixth_Bondary: %f\n'%((mmax-1.)/nscan))
    txt.write('#Seventh_Bondary: %f\n'%((nmax-1.)/nscan))
    #txt.write('#Signal_Model    Significance_after_BDT  Uncertainty_after_BDT   Baseline_Significance   Baseline_Uncertainty    Percentage_improvement\n')
    txt.write('_______________________________________________________________________________________________________________________________________\n\n')

    f_bkg = ROOT.TFile('outputs/model_%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')
    h_bkg = ROOT.TH1F('h_bkg','h_bkg',nscan,0.,1.)
    h_bkg.Sumw2()
    t_bkg.Draw("bdt_score>>h_bkg","weight*1*(235005./769468.)*((%s>=110&&%s<=160)&&!(%s>=120&&%s<=130)&&eventNumber%%2>=0)"%(var,var,var,var))
    nbkg1 ,dbkg1 = hist_integral(h_bkg, 1, hmax-1)
    nbkg2 ,dbkg2 = hist_integral(h_bkg, hmax, imax-1)
    nbkg3 ,dbkg3 = hist_integral(h_bkg, imax, jmax-1)
    nbkg4 ,dbkg4 = hist_integral(h_bkg, jmax, kmax-1)
    nbkg5 ,dbkg5 = hist_integral(h_bkg, kmax, lmax-1)
    nbkg6 ,dbkg6 = hist_integral(h_bkg, lmax, mmax-1)
    nbkg7 ,dbkg7 = hist_integral(h_bkg, mmax, nmax-1)
    nbkg8 ,dbkg8 = hist_integral(h_bkg, nmax, nscan)

    f_sig = ROOT.TFile('outputs/model_%s/sig.root' % (region))
    t_sig = f_sig.Get('test')
    h_sig = ROOT.TH1F('h_sig','h_sig',nscan,0.,1.)
    h_sig.Sumw2()
    t_sig.Draw("bdt_score>>h_sig","weight*1*(%s>=120&&%s<=130&&eventNumber%%2>=0)"%(var,var))
    nsig1, dsig1 = hist_integral(h_sig, 1, hmax-1)
    nsig2, dsig2 = hist_integral(h_sig, hmax, imax-1)
    nsig3, dsig3 = hist_integral(h_sig, imax, jmax-1)
    nsig4, dsig4 = hist_integral(h_sig, jmax, kmax-1)
    nsig5, dsig5 = hist_integral(h_sig, kmax, lmax-1)
    nsig6, dsig6 = hist_integral(h_sig, lmax, mmax-1)
    nsig7, dsig7 = hist_integral(h_sig, mmax, nmax-1)
    nsig8, dsig8 = hist_integral(h_sig, nmax, nscan)

    s1,u1=calc_sig(nsig1,nbkg1,dsig1,dbkg1)
    s2,u2=calc_sig(nsig2,nbkg2,dsig2,dbkg2) 
    s3,u3=calc_sig(nsig3,nbkg3,dsig3,dbkg3)
    s4,u4=calc_sig(nsig4,nbkg4,dsig4,dbkg4)
    s5,u5=calc_sig(nsig5,nbkg5,dsig5,dbkg5)
    s6,u6=calc_sig(nsig6,nbkg6,dsig6,dbkg6)
    s7,u7=calc_sig(nsig7,nbkg7,dsig7,dbkg7)
    s8,u8=calc_sig(nsig8,nbkg8,dsig8,dbkg8)
    s, u = sum_z([s1,s2,s3,s4,s5,s6,s7,s8],[u1,u2,u3,u4,u5,u6,u7,u8])
    #print nsig1, nbkg1, nsig2, nbkg2, nsig3, nbkg3, nsig4, nbkg4, nsig5, nbkg5, nsig6, nbkg6, nsig7, nbkg7, nsig8, nbkg8, s1, s2, s3, s4, s5, s6, s7, s8, abs(u1), abs(u2), abs(u3), abs(u4), abs(u5), abs(u6), abs(u7), abs(u8)
       
    print 'Significance:  %f +- %f'%(s,abs(u))
    txt.write('Significance: %f +- %f\n'%(s,abs(u)))


def categorizing(region,sigs,bkgs,nscan, minN):

    siglist=''
    for sig in sigs:
        if os.path.isfile('outputs/model_%s/%s.root'% (region,sig)): siglist+=' outputs/model_%s/%s.root'% (region,sig)
    os.system("hadd -f outputs/model_%s/sig.root"%(region)+siglist)

    f_sig = ROOT.TFile('outputs/model_%s/sig.root' % (region))
    t_sig = f_sig.Get('test')

    #if not os.path.isfile('outputs/model_%s_%s/bkg.root' % (train_sig,region)):
    bkglist=''
    for bkg in bkgs:
        if os.path.isfile('outputs/model_%s/%s.root'% (region,bkg)): bkglist+=' outputs/model_%s/%s.root'% (region,bkg)
    os.system("hadd -f outputs/model_%s/bkg.root"%(region)+bkglist)
   
    f_bkg = ROOT.TFile('outputs/model_%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig=ROOT.TH1F('h_sig','h_sig',nscan,0,1)
    h_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0,1)

    h_sig.Sumw2()
    h_bkg.Sumw2()


    var='m_mumu'

    t_sig.Draw("bdt_score>>h_sig","weight*1*((%s>=120&&%s<=130)&&(eventNumber%%2>=0))"%(var,var))
    t_bkg.Draw("bdt_score>>h_bkg","weight*1*(235005./769468.)*((%s>=110&&%s<=160)&&!(%s>=120&&%s<=130)&&(eventNumber%%2>=0))"%(var,var,var,var))

    hmax=0
    imax=0
    jmax=0
    kmax=0
    lmax=0
    mmax=0
    nmax=0
    nsig1max=0
    nbkg1max=0
    nsig2max=0
    nbkg2max=0
    nsig3max=0
    nbkg3max=0
    nsig4max=0
    nbkg4max=0
    nsig5max=0
    nbkg5max=0
    nsig6max=0
    nbkg6max=0
    nsig7max=0
    nbkg7max=0
    nsig8max=0
    nbkg8max=0
    s1max=0
    s2max=0
    s3max=0
    s4max=0
    s5max=0
    s6max=0
    s7max=0
    s8max=0
    smax=0

    for k in range(1,nscan+1):

        hmax_low=0
        imax_low=0
        jmax_low=0
        lmax_high=0
        mmax_high=0
        nmax_high=0
        nsig1max_low=0
        nbkg1max_low=0
        nsig2max_low=0
        nbkg2max_low=0
        nsig3max_low=0
        nbkg3max_low=0
        nsig4max_low=0
        nbkg4max_low=0
        nsig5max_high=0
        nbkg5max_high=0
        nsig6max_high=0
        nbkg6max_high=0
        nsig7max_high=0
        nbkg7max_high=0
        nsig8max_high=0
        nbkg8max_high=0
        s1max_low=0
        s2max_low=0
        s3max_low=0
        s4max_low=0
        s5max_high=0
        s6max_high=0
        s7max_high=0
        s8max_high=0        
        smax_low=0
        smax_high=0

        for i in range(1,k+1):

            hmax_low_low=0
            jmax_low_high=0
            nsig1max_low_low=0
            nbkg1max_low_low=0
            nsig2max_low_low=0
            nbkg2max_low_low=0
            nsig3max_low_high=0
            nbkg3max_low_high=0
            nsig4max_low_high=0
            nbkg4max_low_high=0
            s1max_low_low=0
            s2max_low_low=0
            s3max_low_high=0
            s4max_low_high=0
            smax_low_low=0
            smax_low_high=0
                    
            for h in range(1,i+1):

                nsig1, dsig1 = hist_integral(h_sig, 1, h-1)
                nbkg1, dbkg1 = hist_integral(h_bkg, 1, h-1)
                if nbkg1 < minN: continue
                s1, u1 = calc_sig(nsig1, nbkg1, dsig1, dbkg1)
                nsig2, dsig2 = hist_integral(h_sig, h, i-1)
                nbkg2, dbkg2 = hist_integral(h_bkg, h, i-1)
                if nbkg2 < minN: continue
                s2, u2 = calc_sig(nsig2, nbkg2, dsig2, dbkg2)
                s_low_low, u_low_low = sum_z([s1,s2],[u1,u2])
                if s_low_low > smax_low_low:
                    smax_low_low = s_low_low
                    hmax_low_low = h
                    nsig1max_low_low = nsig1
                    nbkg1max_low_low = nbkg1
                    nsig2max_low_low = nsig2
                    nbkg2max_low_low = nbkg2
                    s1max_low_low = s1
                    s2max_low_low = s2

            for j in range(i,k+1):

                nsig3, dsig3 = hist_integral(h_sig, i, j-1)
                nbkg3, dbkg3 = hist_integral(h_bkg, i, j-1)
                if nbkg3 < minN: continue
                s3, u3 = calc_sig(nsig3, nbkg3, dsig3, dbkg3)
                nsig4, dsig4 = hist_integral(h_sig, j, k-1)
                nbkg4, dbkg4 = hist_integral(h_bkg, j, k-1)
                if nbkg4 < minN: continue
                s4, u4 = calc_sig(nsig4, nbkg4, dsig4, dbkg4)
                s_low_high, u_low_high = sum_z([s3,s4],[u3,u4])
                if s_low_high > smax_low_high:
                    smax_low_high = s_low_high
                    jmax_low_high = j
                    nsig3max_low_high = nsig3
                    nbkg3max_low_high = nbkg3
                    nsig4max_low_high = nsig4
                    nbkg4max_low_high = nbkg4
                    s3max_low_high = s3
                    s4max_low_high = s4

            s_low = sqrt(smax_low_low**2 + smax_low_high**2)
            if s_low>smax_low:
                smax_low=s_low
                jmax_low=jmax_low_high
                imax_low=i
                hmax_low=hmax_low_low
                nsig4max_low=nsig4max_low_high
                nbkg4max_low=nbkg4max_low_high
                nsig3max_low=nsig3max_low_high
                nbkg3max_low=nbkg3max_low_high
                nsig2max_low=nsig2max_low_low
                nbkg2max_low=nbkg2max_low_low
                nsig1max_low=nsig1max_low_low
                nbkg1max_low=nbkg1max_low_low
                s4max_low=s4max_low_high
                s3max_low=s3max_low_high
                s2max_low=s2max_low_low
                s1max_low=s1max_low_low

        for m in range(k,nscan+1):

            lmax_high_low=0
            nmax_high_high=0
            nsig5max_high_low=0
            nbkg5max_high_low=0
            nsig6max_high_low=0
            nbkg6max_high_low=0
            nsig7max_high_high=0
            nbkg7max_high_high=0
            nsig8max_high_high=0
            nbkg8max_high_high=0
            s5max_high_low=0
            s6max_high_low=0
            s7max_high_high=0
            s8max_high_high=0
            smax_high_low=0
            smax_high_high=0

            for l in range(k,m+1):

                nsig5, dsig5 = hist_integral(h_sig, k, l-1)
                nbkg5, dbkg5 = hist_integral(h_bkg, k, l-1)
                if nbkg5 < minN: continue
                s5, u5 = calc_sig(nsig5, nbkg5, dsig5, dbkg5)
                nsig6, dsig6 = hist_integral(h_sig, l, m-1)
                nbkg6, dbkg6 = hist_integral(h_bkg, l, m-1)
                if nbkg6 < minN: continue
                s6, u6 = calc_sig(nsig6, nbkg6, dsig6, dbkg6)
                s_high_low, u_high_low = sum_z([s5,s6],[u5,u6])
                if s_high_low > smax_high_low:
                    smax_high_low = s_high_low
                    lmax_high_low = l
                    nsig5max_high_low = nsig5
                    nbkg5max_high_low = nbkg5
                    nsig6max_high_low = nsig6
                    nbkg6max_high_low = nbkg6
                    s5max_high_low = s5
                    s6max_high_low = s6

            for n in range(m,nscan+1):

                nsig7, dsig7 = hist_integral(h_sig, m, n-1)
                nbkg7, dbkg7 = hist_integral(h_bkg, m, n-1)
                if nbkg7 < minN: continue
                s7, u7 = calc_sig(nsig7, nbkg7, dsig7, dbkg7)
                nsig8, dsig8 = hist_integral(h_sig, n, nscan)
                nbkg8, dbkg8 = hist_integral(h_bkg, n, nscan)
                if nbkg8 < minN: continue
                s8, u8 = calc_sig(nsig8, nbkg8, dsig8, dbkg8)
                s_high_high, u_high_high = sum_z([s7,s8],[u7,u8])
                if s_high_high > smax_high_high:
                    smax_high_high = s_high_high
                    nmax_high_high = n
                    nsig7max_high_high = nsig7
                    nbkg7max_high_high = nbkg7
                    nsig8max_high_high = nsig8
                    nbkg8max_high_high = nbkg8
                    s7max_high_high = s7
                    s8max_high_high = s8

            s_high = sqrt(smax_high_low**2 + smax_high_high**2)
            if s_high>smax_high:
                smax_high=s_high
                nmax_high=nmax_high_high
                mmax_high=m
                lmax_high=lmax_high_low
                nsig5max_high=nsig5max_high_low
                nbkg5max_high=nbkg5max_high_low
                nsig6max_high=nsig6max_high_low
                nbkg6max_high=nbkg6max_high_low
                nsig7max_high=nsig7max_high_high
                nbkg7max_high=nbkg7max_high_high
                nsig8max_high=nsig8max_high_high
                nbkg8max_high=nbkg8max_high_high
                s5max_high=s5max_high_low
                s6max_high=s6max_high_low
                s7max_high=s7max_high_high
                s8max_high=s8max_high_high

        s=sqrt(smax_low**2+smax_high**2)
        if s>smax:
            smax=s
            hmax=hmax_low
            imax=imax_low
            jmax=jmax_low
            kmax=k
            lmax=lmax_high
            mmax=mmax_high
            nmax=nmax_high
            nsig1max=nsig1max_low
            nbkg1max=nbkg1max_low
            nsig2max=nsig2max_low
            nbkg2max=nbkg2max_low
            nsig3max=nsig3max_low
            nbkg3max=nbkg3max_low
            nsig4max=nsig4max_low
            nbkg4max=nbkg4max_low
            nsig5max=nsig5max_high
            nbkg5max=nbkg5max_high
            nsig6max=nsig6max_high
            nbkg6max=nbkg6max_high
            nsig7max=nsig7max_high
            nbkg7max=nbkg7max_high
            nsig8max=nsig8max_high
            nbkg8max=nbkg8max_high
            s1=s1max_low
            s2=s2max_low
            s3=s3max_low
            s4=s4max_low
            s5=s5max_high
            s6=s6max_high
            s7=s7max_high
            s8=s8max_high
            

    #print nsig1max, nbkg1max, nsig2max, nbkg2max, nsig3max, nbkg3max, nsig4max, nbkg4max, nsig5max, nbkg5max, nsig6max, nbkg6max, nsig7max, nbkg7max, nsig8max, nbkg8max, s1max, s2max, s3max, s4max, s5max, s6max, s7max, s8max

    print '========================================================================='
    print 'The maximal significance:  %f' %(smax)
    print 'First boundary:  %f,    Second boundary  %f,    Third boundary  %f,    Fourth boundary  %f,    Fifth boundary  %f,    Sixth boundary  %f,    Seventh boundary  %f' %((hmax-1.)/nscan, (imax-1.)/nscan, (jmax-1.)/nscan, (kmax-1.)/(nscan), (lmax-1.)/(nscan), (mmax-1.)/(nscan), (nmax-1.)/(nscan))
    print '========================================================================='

    return hmax, imax, jmax, kmax, lmax, mmax, nmax


def main():

    ROOT.gROOT.SetBatch(True)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data']




    if args.region == 'two_jet':
        region = 'two_jet'
    elif args.region == 'one_jet':
        region = 'one_jet'
    else:
        region = 'zero_jet'

    nscan=args.nscan

    CrossSections={}

    hmax, imax, jmax, kmax, lmax, mmax, nmax=categorizing(region,sigs,bkgs, nscan, args.minN)

    gettingsig(region,sigs,bkgs, hmax,imax,jmax,kmax,lmax,mmax,nmax,nscan)

    return



if __name__ == '__main__':
    main()
