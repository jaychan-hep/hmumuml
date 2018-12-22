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
    return  parser.parse_args()

def calc_sig(sig,bkg,d_sig,d_bkg):

  ntot = sig + bkg

  if(sig <= 0): return 0, 0
  if(bkg <= 0): return 0, 0

  #Counting experiment
  significance = sqrt(2*((ntot*log(ntot/bkg)) - sig))

  #error on significance
  numer = log(1+(sig/bkg)) - (sig/bkg)
  #delta_b = sig/sqrt(sig/0.00001)
  #if (isTest) delta_b = n_data/sqrt(n_data/(4./45.))
  uncert = (numer/significance) * d_bkg

  #safety feature to prevent low number of events
  #if(n_data < 0.8) {significance = 0; uncert=0;}
  return significance, uncert

def hist_integral(hist,i,j):
    if i>j: n=0
    else: n = hist.Integral(i,j)
    return n

def gettingsig(region,sigs,bkgs, hmax,imax,jmax,kmax,lmax,mmax,nmax,nscan):

    print 'Evaluating significances for all of the signal models...'

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
    w_bkg=5
    d_bkg1=ROOT.Double(0)
    if (1>hmax-1): nbkg1=0
    else: nbkg1=h_bkg.IntegralAndError(1,hmax-1,d_bkg1)
    d_bkg2=ROOT.Double(0)
    if (hmax>imax-1): nbkg2=0
    else: nbkg2=h_bkg.IntegralAndError(hmax,imax-1,d_bkg2)
    d_bkg3=ROOT.Double(0)
    if (imax>jmax-1): nbkg3=0
    else: nbkg3=h_bkg.IntegralAndError(imax,jmax-1,d_bkg3)
    d_bkg4=ROOT.Double(0)
    if (jmax>kmax-1): nbkg4=0
    else: nbkg4=h_bkg.IntegralAndError(jmax,kmax-1,d_bkg4)
    d_bkg5=ROOT.Double(0)
    if (kmax>lmax-1): nbkg5=0
    else: nbkg5=h_bkg.IntegralAndError(kmax,lmax-1,d_bkg5)
    d_bkg6=ROOT.Double(0)
    if (lmax>mmax-1): nbkg6=0
    else: nbkg6=h_bkg.IntegralAndError(lmax,mmax-1,d_bkg6)
    d_bkg7=ROOT.Double(0)
    if (mmax>nmax-1): nbkg7=0
    else: nbkg7=h_bkg.IntegralAndError(mmax,nmax-1,d_bkg7)
    d_bkg8=ROOT.Double(0)
    if (nmax>nscan): nbkg8=0
    else: nbkg8=h_bkg.IntegralAndError(nmax,nscan,d_bkg8)
    #nbkg=baseline_yield['bkg_tot']

    f_sig = ROOT.TFile('outputs/model_%s/sig.root' % (region))
    t_sig = f_sig.Get('test')
    h_sig = ROOT.TH1F('h_sig','h_sig',nscan,0.,1.)
    h_sig.Sumw2()
    t_sig.Draw("bdt_score>>h_sig","weight*1*(%s>=120&&%s<=130&&eventNumber%%2>=0)"%(var,var))
    nsig1=hist_integral(h_sig,1,hmax-1)
    nsig2=hist_integral(h_sig,hmax,imax-1)
    nsig3=hist_integral(h_sig,imax,jmax-1)
    nsig4=hist_integral(h_sig,jmax,kmax-1)
    nsig5=hist_integral(h_sig,kmax,lmax-1)
    nsig6=hist_integral(h_sig,lmax,mmax-1)
    nsig7=hist_integral(h_sig,mmax,nmax-1)
    nsig8=hist_integral(h_sig,nmax,nscan)

    #nsig1=nsig1*443.0/435.7
    #nsig2=nsig2*443.0/435.7
    #nsig3=nsig3*443.0/435.7
    #nsig4=nsig4*443.0/435.7
    #nsig5=nsig5*443.0/435.7

    #nbkg1=nbkg1*764530./769468.
    #nbkg2=nbkg2*764530./769468.
    #nbkg3=nbkg3*764530./769468.
    #nbkg4=nbkg4*764530./769468.
    #nbkg5=nbkg5*764530./769468.

    s1,u1=calc_sig(nsig1,nbkg1,0,d_bkg1)
    s2,u2=calc_sig(nsig2,nbkg2,0,d_bkg2) 
    s3,u3=calc_sig(nsig3,nbkg3,0,d_bkg3)
    s4,u4=calc_sig(nsig4,nbkg4,0,d_bkg4)
    s5,u5=calc_sig(nsig5,nbkg5,0,d_bkg5)
    s6,u6=calc_sig(nsig6,nbkg6,0,d_bkg6)
    s7,u7=calc_sig(nsig7,nbkg7,0,d_bkg7)
    s8,u8=calc_sig(nsig8,nbkg8,0,d_bkg8)
    s=sqrt(s1**2+s2**2+s3**2+s4**2+s5**2+s6**2+s7**2+s8**2)
    u=sqrt((s1*u1)**2+(s2*u2)**2+(s3*u3)**2+(s4*u4)**2+(s5*u5)**2+(s6*u6)**2+(s7*u7)**2+(s8*u8)**2)/s
    #nsig=baseline_yield[sig]*CrossSections[sig.replace('zp2hdmbbmzp','').replace('mA',',')]
    #s0,u0=calc_sig(nsig,nbkg,1,1)
    print nsig1, nbkg1, nsig2, nbkg2, nsig3, nbkg3, nsig4, nbkg4, nsig5, nbkg5, nsig6, nbkg6, nsig7, nbkg7, nsig8, nbkg8, s1, s2, s3, s4, s5, s6, s7, s8, abs(u1), abs(u2), abs(u3), abs(u4), abs(u5), abs(u6), abs(u7), abs(u8)
       
    #print 'Significance of %s:  %f +- %f         Baseline significance:  %f +- %f        Imporvement: %f %%'%(sig,s,abs(u),s0,abs(u0),(s-s0)/s0*100 if s0 !=0 else 0)
    #txt.write('%s  %f  %f  %f  %f  %f\n'%(sig,s,abs(u),s0,abs(u0),(s-s0)/s0*100 if s0 !=0 else 0))
    print 'Significance:  %f +- %f'%(s,abs(u))
    txt.write('%f  %f\n'%(s,abs(u)))


def categorizing(region,sigs,bkgs,nscan):

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

        for j in range(1,k+1):
            for i in range(1,j+1):
                for h in range(1,i+1):
                    nsig4=hist_integral(h_sig,j,k-1)
                    d_bkg4=ROOT.Double(0)
                    if (j>k-1): nbkg4=0
                    else: nbkg4=h_bkg.IntegralAndError(j,k-1,d_bkg4)
                    s4,u4=calc_sig(nsig4,nbkg4,0,d_bkg4)
                    if nbkg4<2: s4,u4=0,0
                    nsig3=hist_integral(h_sig,i,j-1)
                    d_bkg3=ROOT.Double(0)
                    if (i>j-1): nbkg3=0
                    else: nbkg3=h_bkg.IntegralAndError(i,j-1,d_bkg3)
                    s3,u3=calc_sig(nsig3,nbkg3,0,d_bkg3)
                    if nbkg3<2: s3,u3=0,0
                    nsig2=hist_integral(h_sig,h,i-1)
                    d_bkg2=ROOT.Double(0)
                    if (h>i-1): nbkg2=0
                    else: nbkg2=h_bkg.IntegralAndError(h,i-1,d_bkg2)
                    s2,u2=calc_sig(nsig2,nbkg2,0,d_bkg2)
                    if nbkg2<2: s2,u2=0,0
                    nsig1=hist_integral(h_sig,1,h-1)
                    d_bkg1=ROOT.Double(0)
                    if (1>h-1): nbkg1=0
                    else: nbkg1=h_bkg.IntegralAndError(1,h-1,d_bkg1)
                    s1,u1=calc_sig(nsig1,nbkg1,0,d_bkg1)
                    if nbkg1<2: s1,u1=0,0
                    s_low=sqrt(s1**2+s2**2+s3**2+s4**2)
                    if s_low>smax_low: 
                        smax_low=s_low
                        jmax_low=j
                        imax_low=i
                        hmax_low=h
                        nsig4max_low=nsig4
                        nbkg4max_low=nbkg4
                        nsig3max_low=nsig3
                        nbkg3max_low=nbkg3
                        nsig2max_low=nsig2
                        nbkg2max_low=nbkg2
                        nsig1max_low=nsig1
                        nbkg1max_low=nbkg1
                        s4max_low=s4
                        s3max_low=s3
                        s2max_low=s2
                        s1max_low=s1

        for l in range(k,nscan+1):
            for m in range(l,nscan+1):
                for n in range(m,nscan+1):
                    nsig5=hist_integral(h_sig,k,l-1)
                    d_bkg5=ROOT.Double(0)
                    if (k>l-1): nbkg5=0
                    else: nbkg5=h_bkg.IntegralAndError(k,l-1,d_bkg5)
                    s5,u5=calc_sig(nsig5,nbkg5,0,d_bkg5)
                    if nbkg5<2: s5,u5=0,0
                    nsig6=hist_integral(h_sig,l,m-1)
                    d_bkg6=ROOT.Double(0)
                    if (l>m-1): nbkg6=0
                    else: nbkg6=h_bkg.IntegralAndError(l,m-1,d_bkg6)
                    s6,u6=calc_sig(nsig6,nbkg6,0,d_bkg6)
                    if nbkg6<2: s6,u6=0,0
                    nsig7=hist_integral(h_sig,m,n-1)
                    d_bkg7=ROOT.Double(0)
                    if (m>n-1): nbkg7=0
                    else: nbkg7=h_bkg.IntegralAndError(m,n-1,d_bkg7)
                    s7,u7=calc_sig(nsig7,nbkg7,0,d_bkg7)
                    if nbkg7<2: s7,u7=0,0
                    nsig8=hist_integral(h_sig,n,nscan)
                    d_bkg8=ROOT.Double(0)
                    if (n>nscan): nbkg8=0
                    else: nbkg8=h_bkg.IntegralAndError(n,nscan,d_bkg8)
                    s8,u8=calc_sig(nsig8,nbkg8,0,d_bkg8)
                    if nbkg8<2: s8,u8=0,0
                    s_high=sqrt(s5**2+s6**2+s7**2+s8**2)
                    if s_high>smax_high: 
                        smax_high=s_high
                        lmax_high=l
                        mmax_high=m
                        nmax_high=n
                        nsig5max_high=nsig5
                        nbkg5max_high=nbkg5
                        nsig6max_high=nsig6
                        nbkg6max_high=nbkg6
                        nsig7max_high=nsig7
                        nbkg7max_high=nbkg7
                        nsig8max_high=nsig8
                        nbkg8max_high=nbkg8
                        s5max_high=s5
                        s6max_high=s6
                        s7max_high=s7
                        s8max_high=s8

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
            

    print nsig1max, nbkg1max, nsig2max, nbkg2max, nsig3max, nbkg3max, nsig4max, nbkg4max, nsig5max, nbkg5max, nsig6max, nbkg6max, nsig7max, nbkg7max, nsig8max, nbkg8max, s1max, s2max, s3max, s4max, s5max, s6max, s7max, s8max

    print '========================================================================='
    print 'The maximal significance:  %f' %(smax)
    print 'First boundary:  %f,    Second boundary  %f,    Third boundary  %f,    Fourth boundary  %f,    Fifth boundary  %f,    Sixth boundary  %f,    Seventh boundary  %f' %((hmax-1.)/nscan,(imax-1.)/nscan,(jmax-1.)/nscan,(kmax-1.)/(nscan),(lmax-1.)/(nscan),(mmax-1.)/(nscan),(nmax-1.)/(nscan))
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

    hmax, imax, jmax, kmax, lmax, mmax, nmax=categorizing(region,sigs,bkgs, nscan)

    gettingsig(region,sigs,bkgs, hmax,imax,jmax,kmax,lmax,mmax,nmax,nscan)

    return



if __name__ == '__main__':
    main()
