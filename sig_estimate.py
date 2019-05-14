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
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'all_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-n','--nscan',type=int,default=2000,help='number of scan.')
    parser.add_argument('--minN',type=float,default=50,help='minimum number of events in mass window')
    return  parser.parse_args()

def calc_sig(sig, bkg,s_err,b_err):

  ntot = sig + bkg

  if(sig <= 0): return 0, 0
  if(bkg <= 0): return 0, 0

  #Counting experiment
  #significance = sqrt(2*((ntot*log(ntot/bkg)) - sig))
  significance = sig / sqrt(bkg)

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

def categorizing(region,sigs,bkgs,nscan, minN):

    #siglist=''
    #for sig in sigs:
    #    if os.path.isfile('outputs/model_%s/%s.root'% (region,sig)): siglist+=' outputs/model_%s/%s.root'% (region,sig)
    #os.system("hadd -f outputs/model_%s/sig.root"%(region)+siglist)

    f_sig = ROOT.TFile('outputs/model_%s/sig.root' % (region))
    t_sig = f_sig.Get('test')

    #if not os.path.isfile('outputs/model_%s_%s/bkg.root' % (train_sig,region)):
    #bkglist=''
    #for bkg in bkgs:
    #    if os.path.isfile('outputs/model_%s/%s.root'% (region,bkg)): bkglist+=' outputs/model_%s/%s.root'% (region,bkg)
    #os.system("hadd -f outputs/model_%s/bkg.root"%(region)+bkglist)
   
    f_bkg = ROOT.TFile('outputs/model_%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_sig=ROOT.TH1F('h_sig','h_sig',nscan,0,1)
    h_bkg=ROOT.TH1F('h_bkg','h_bkg',nscan,0,1)

    h_sig.Sumw2()
    h_bkg.Sumw2()


    var='m_mumu'

    t_sig.Draw("bdt_score_t>>h_sig","weight*1*((%s>=120&&%s<=130)&&(eventNumber%%2>=0))"%(var,var))
    t_bkg.Draw("bdt_score_t>>h_bkg","weight*1*(235005./769468.)*((%s>=110&&%s<=160)&&!(%s>=120&&%s<=130)&&(eventNumber%%2>=0))"%(var,var,var,var))


    s2 = 0
    nsig = 0
    nbkg = 0
    ncat = 0

    for i in range(nscan):

        nsig += hist_integral(h_sig, nscan-i, nscan-i)[0]
        nbkg += hist_integral(h_bkg, nscan-i, nscan-i)[0]
        if nbkg > minN:
            s2 += nsig**2/nbkg
            nsig = 0
            nbkg = 0
            ncat += 1
            print 'INFO: Found one boundary: %f' %((nscan-i-1.)/nscan)

    #print nsig1max, nbkg1max, nsig2max, nbkg2max, nsig3max, nbkg3max, nsig4max, nbkg4max, nsig5max, nbkg5max, nsig6max, nbkg6max, nsig7max, nbkg7max, nsig8max, nbkg8max, s1max, s2max, s3max, s4max, s5max, s6max, s7max, s8max

    print '========================================================================='
    print 'The estimate significance:  %f' %(sqrt(s2))
    print 'Number of categories:  %f' %(ncat)
    print '========================================================================='

    return


def main():

    ROOT.gROOT.SetBatch(True)

    args=getArgs()

    sigs = ['ggF','VBF','VH','ttH']

    bkgs = ['data']



    region = args.region

    nscan=args.nscan

    CrossSections={}

    categorizing(region,sigs,bkgs, nscan, args.minN)

    return



if __name__ == '__main__':
    main()
