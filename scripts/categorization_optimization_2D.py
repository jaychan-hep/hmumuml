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
from math import sqrt, log
import math
import sys
import time
import json

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

#def hist_integral(hist,i,j):
#    err = Double()
#    if i>j:
#        n=0
#        err=0
#    else: n = hist.IntegralAndError(i,j,err)
#    return n, err

def hist_integral(hist,i,j,k=-999,l=-999):
    err = Double()
    if i>j or k>l:
        n=0
        err=0
    elif k == -999 and l == -999: n = hist.IntegralAndError(i,j,err)
    else: n = hist.IntegralAndError(i,j,k,l,err)
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

def fit_BDT(hist, fitboundary, fitbin, function, fillbin, sf):

    n_events = hist.Integral()

    bdt = RooRealVar("bdt","bdt",fitboundary,1)

    #lam = RooRealVar("lam", "lam", 10, -100, 20)
    a = RooRealVar("a", "a", -100, 100)
    b = RooRealVar("b", "b", -100, 100)
    c = RooRealVar("c", "c", -100, 100)
    d = RooRealVar("d", "d", -100, 100)
    e = RooRealVar("e", "e", -100, 100)
    f = RooRealVar("f", "f", -100, 100)

    pdf = {}
    pdf['power'] = RooGenericPdf("pdf","PDF for BDT distribution", "@0**@1", RooArgList(bdt, a)) #power

    # EPoly Family
    pdf['Exp'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0)", RooArgList(bdt, a))
    pdf['Epoly2'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2)", RooArgList(bdt, a, b))
    pdf['Epoly3'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3)", RooArgList(bdt, a, b, c))
    pdf['Epoly4'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4)", RooArgList(bdt, a, b, c, d))
    pdf['Epoly5'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5)", RooArgList(bdt, a, b, c, d, e))
    pdf['Epoly6'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5+@6*@0**6)", RooArgList(bdt, a, b, c, d, e, f))

    #pdf['Poly6'] = RooGenericPdf("pdf","PDF for BDT distribution", "1+@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5+@6*@0**6", RooArgList(bdt, a, b, c, d, e, f))
    #pdf['ExPoly6'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(-@0*@1)*(1+@2*@0+@3*@0**2+@4*@0**3+@5*@0**4+@6*@0**5+@7*@0**6)", RooArgList(bdt, lam, a, b, c, d, e, f))

    data = RooDataHist("dh", "dh", RooArgList(bdt), hist)
    #data.Print("all")


    pdf[function].fitTo(data, RooFit.Verbose(False), RooFit.PrintLevel(-1))

    #lam.Print()
    a.Print()
    b.Print()
    c.Print()
    d.Print()
    e.Print()
    f.Print()

    frame = bdt.frame()
    data.plotOn(frame)
    pdf[function].plotOn(frame)

    frame.Draw()

    dof = {'power': 2, 'Exp': 2, 'Epoly2': 3, 'Epoly3': 4, 'Epoly4': 5, 'Epoly5': 6, 'Epoly6': 7}
    reduced_chi_square = frame.chiSquare(dof[function])
    probability = TMath.Prob(frame.chiSquare(dof[function]) * (fitbin - dof[function]), fitbin - dof[function])
    print 'chi square:', reduced_chi_square
    print 'probability: ', probability

    #raw_input('Press enter to continue')

    # fill the fitted pdf into a histogram
    hfit = TH1F('hfit', 'hfit', fillbin, fitboundary, 1.)
    pdf[function].fillHistogram(hfit, RooArgList(bdt), n_events*sf)

    return hfit


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


    # prepare some settings for bkg histograms
    fitboundary = 0.5
    nbin_right = int(nscanvbf*(1-fitboundary))
    nbin_left = int(nscanvbf*fitboundary)
    fitbin = nbin_right  #100

    h_bkg_raw = TH2F('h_bkg_raw', 'h_bkg_raw', nscanvbf, 0., 1., nscan, 0., 1.)
    h_bkg_right = TH1F('h_bkg_right', 'h_bkg_right', fitbin, fitboundary, 1.)
    h_bkg_left = TH1F('h_bkg_left', 'h_bkg_left', nbin_left, 0., fitboundary)
    h_bkg_raw.Sumw2()
    h_bkg_right.Sumw2()
    h_bkg_left.Sumw2()

    # filling bkg histograms
    t_bkg.Draw("bdt_score%s:bdt_score_VBF%s>>h_bkg_raw"%('_t' if transform else '','_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1, n_fold, fold if n_fold != 1 else 1))
    t_bkg.Draw("bdt_score_VBF%s>>h_bkg_right"%('_t' if transform else ''),"weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d))"%(n_fold,fold if n_fold != 1 else 1))
    t_bkg.Draw("bdt_score_VBF%s>>h_bkg_left"%('_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(bdt_score_VBF_t<%f&&bdt_score_VBF_t>=0&&(eventNumber%%%d!=%d)))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,fitboundary, n_fold,fold if n_fold != 1 else 1))

    sf = (n_fold/(n_fold-1.) if n_fold != 1 else 1)*0.2723

    # fit a function to the bkg BDT distribution!
    h_bkg_right = fit_BDT(h_bkg_right, fitboundary, fitbin, 'Epoly2', nbin_right, sf)
    #h_bkg_right.Scale(sf) # uncomment this line (and comment out the line above) to use the original histogram instead of the fitted histogram


    # merge the non-fitted histogram with the fitted histogram
    h_bkg = TH1F('h_bkg','h_bkg',nscanvbf,0.,1.)
    h_bkg.Sumw2()
    h_bkg_list = TList()
    h_bkg_list.Add(h_bkg_left)
    h_bkg_list.Add(h_bkg_right)
    h_bkg.Merge(h_bkg_list)

    #================================================
    # categorization for VBF categories. Will scan all of the possibilities of the BDT boundaries and get the one that gives the largest significance

    jv1max, jv2max, jv3max, jv4max, jv5max, jv6max, jv7max, jv8max, jv9max, jv10max, jv11max, jv12max, jv13max, jv14max, jv15max, svmax = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    print 'INFO: scanning all of the possibilities for VBF categories...'
    start_time = time.time()
    jv8stop = 0
    for jv8 in range(vb, nscanvbf+1 if nvbf >= 9 else vb+1):

        elapsed_time = time.time() - start_time
        tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print '%s |%s%s| %d/100'%(tt, '>'*(jv8*100/nscanvbf), ' '*(100-jv8*100/nscanvbf), jv8*100/nscanvbf)
        sys.stdout.write("\033[F") # Cursor up one line

        jv1max_, jv2max_, jv3max_, jv4max_, jv5max_, jv6max_, jv7max_, svmax_1_8 = 0, 0, 0, 0, 0, 0, 0, 0
        jv4stop = 0
        for jv4 in range(vb, jv8+1 if nvbf >= 13 else vb+1):

            jv1max__, jv2max__, jv3max__, jv5max__, jv6max__, jv7max__, svmax_1_4, svmax_5_8 = 0, 0, 0, 0, 0, 0, 0, 0
            jv2stop = 0
            for jv2 in range(vb, jv8+1 if nvbf >= 15 else vb+1):

                jv1max___, jv3max___, svmax_1_2, svmax_3_4 = 0, 0, 0, 0
                jv1stop = 0
                for jv1 in range(vb, jv2+1 if nvbf >= 16 else vb+1):
                    nsigv1, dsigv1 = hist_integral(h_sig, vb, jv1-1)
                    nbkgv1, dbkgv1 = hist_integral(h_bkg, vb, jv1-1)
                    if nvbf >= 16 and nbkgv1 < minN: continue
                    sv1, uv1 = calc_sig(nsigv1, nbkgv1, dsigv1, dbkgv1)
                    nsigv2, dsigv2 = hist_integral(h_sig, jv1, jv2-1)
                    nbkgv2, dbkgv2 = hist_integral(h_bkg, jv1, jv2-1)
                    if nvbf >= 15 and nbkgv2 < minN: break
                    sv2, uv2 = calc_sig(nsigv2, nbkgv2, dsigv2, dbkgv2)
                    sv_1_2, _ = sum_z([[sv1,uv1],[sv2,uv2]])
                    if sv_1_2 >= svmax_1_2:
                        jv1stop = 0
                        svmax_1_2 = sv_1_2
                        jv1max___ = jv1
                    else:
                        jv1stop += 1
                    if jv1stop == earlystop: break

                if jv1max___ == 0: continue

                jv3stop = 0
                for jv3 in range(jv2, jv4+1 if nvbf >= 14 else vb+1):
                    nsigv3, dsigv3 = hist_integral(h_sig, jv2, jv3-1)
                    nbkgv3, dbkgv3 = hist_integral(h_bkg, jv2, jv3-1)
                    if nvbf >= 14 and nbkgv3 < minN: continue
                    sv3, uv3 = calc_sig(nsigv3, nbkgv3, dsigv3, dbkgv3)
                    nsigv4, dsigv4 = hist_integral(h_sig, jv3, jv4-1)
                    nbkgv4, dbkgv4 = hist_integral(h_bkg, jv3, jv4-1)
                    if nvbf >= 13 and nbkgv4 < minN: break
                    sv4, uv4 = calc_sig(nsigv4, nbkgv4, dsigv4, dbkgv4)
                    sv_3_4, _ = sum_z([[sv3,uv3],[sv4,uv4]])
                    if sv_3_4 >= svmax_3_4:
                        jv3stop = 0
                        svmax_3_4 = sv_3_4
                        jv3max___ = jv3
                    else:
                        jv3stop += 1
                    if jv3stop == earlystop: break

                if jv3max___ == 0: break

                sv_1_4 = sqrt(svmax_1_2**2 + svmax_3_4**2)
                if sv_1_4 >= svmax_1_4:
                    jv2stop = 0
                    svmax_1_4 = sv_1_4
                    jv1max__, jv2max__, jv3max__ = jv1max___, jv2, jv3max___
                else:
                    jv2stop += 1
                if jv2stop == earlystop: break

            if jv2max__ == 0: continue

            jv6stop = 0
            for jv6 in range(jv4, jv8+1 if nvbf >= 11 else vb+1):

                jv5max___, jv7max___, svmax_5_6, svmax_7_8 = 0, 0, 0, 0
                jv5stop = 0
                for jv5 in range(jv4, jv6+1 if nvbf >= 12 else vb+1):
                    nsigv5, dsigv5 = hist_integral(h_sig, jv4, jv5-1)
                    nbkgv5, dbkgv5 = hist_integral(h_bkg, jv4, jv5-1)
                    if nvbf >= 12 and nbkgv5 < minN: continue
                    sv5, uv5 = calc_sig(nsigv5, nbkgv5, dsigv5, dbkgv5)
                    nsigv6, dsigv6 = hist_integral(h_sig, jv5, jv6-1)
                    nbkgv6, dbkgv6 = hist_integral(h_bkg, jv5, jv6-1)
                    if nvbf >= 11 and nbkgv6 < minN: break
                    sv6, uv6 = calc_sig(nsigv6, nbkgv6, dsigv6, dbkgv6)
                    sv_5_6, _ = sum_z([[sv5,uv5],[sv6,uv6]])
                    if sv_5_6 >= svmax_5_6:
                        jv5stop = 0
                        svmax_5_6 = sv_5_6
                        jv5max___ = jv5
                    else:
                        jv5stop += 1
                    if jv5stop == earlystop: break

                if jv5max___ == 0: continue

                jv7stop = 0
                for jv7 in range(jv6, jv8+1 if nvbf >= 10 else vb+1):
                    nsigv7, dsigv7 = hist_integral(h_sig, jv6, jv7-1)
                    nbkgv7, dbkgv7 = hist_integral(h_bkg, jv6, jv7-1)
                    if nvbf >= 10 and nbkgv7 < minN: continue
                    sv7, uv7 = calc_sig(nsigv7, nbkgv7, dsigv7, dbkgv7)
                    nsigv8, dsigv8 = hist_integral(h_sig, jv7, jv8-1)
                    nbkgv8, dbkgv8 = hist_integral(h_bkg, jv7, jv8-1)
                    if nvbf >= 9 and nbkgv8 < minN: break
                    sv8, uv8 = calc_sig(nsigv8, nbkgv8, dsigv8, dbkgv8)
                    sv_7_8, _ = sum_z([[sv7,uv7],[sv8,uv8]])
                    if sv_7_8 >= svmax_7_8:
                        jv7stop = 0
                        svmax_7_8 = sv_7_8
                        jv7max___ = jv7
                    else:
                        jv7stop += 1
                    if jv7stop == earlystop: break

                if jv7max___ == 0: break

                sv_5_8 = sqrt(svmax_5_6**2 + svmax_7_8**2)
                if sv_5_8 >= svmax_5_8:
                    jv6stop = 0
                    svmax_5_8 = sv_5_8
                    jv5max__, jv6max__, jv7max__ = jv5max___, jv6, jv7max___
                else:
                    jv6stop += 1
                if jv6stop == earlystop: break

            if jv6max__ == 0: break

            sv_1_8 = sqrt(svmax_1_4**2 + svmax_5_8**2)
            if sv_1_8 >= svmax_1_8:
                jv4stop = 0
                svmax_1_8 = sv_1_8
                jv1max_, jv2max_, jv3max_, jv4max_, jv5max_, jv6max_, jv7max_ = jv1max__, jv2max__, jv3max__, jv4, jv5max__, jv6max__, jv7max__
            else:
                jv4stop += 1
            if jv4stop == earlystop: break

        if jv4max_ == 0: continue
        ################################

        jv9max_, jv10max_, jv11max_, jv12max_, jv13max_, jv14max_, jv15max_, svmax_9_16 = 0, 0, 0, 0, 0, 0, 0, 0
        jv12stop = 0
        for jv12 in range(jv8, nscanvbf+1 if nvbf >= 5 else vb+1):

            jv9max__, jv10max__, jv11max__, jv13max__, jv14max__, jv15max__, svmax_9_12, svmax_13_16 = 0, 0, 0, 0, 0, 0, 0, 0
            jv10stop = 0
            for jv10 in range(jv8, jv12+1 if nvbf >= 7 else vb+1):

                jv9max___, jv11max___, svmax_9_10, svmax_11_12 = 0, 0, 0, 0
                jv9stop = 0
                for jv9 in range(jv8, jv10+1 if nvbf >= 8 else vb+1):
                    nsigv9, dsigv9 = hist_integral(h_sig, jv8, jv9-1)
                    nbkgv9, dbkgv9 = hist_integral(h_bkg, jv8, jv9-1)
                    if nvbf >= 8 and nbkgv9 < minN: continue
                    sv9, uv9 = calc_sig(nsigv9, nbkgv9, dsigv9, dbkgv9)
                    nsigv10, dsigv10 = hist_integral(h_sig, jv9, jv10-1)
                    nbkgv10, dbkgv10 = hist_integral(h_bkg, jv9, jv10-1)
                    if nvbf >= 7 and nbkgv10 < minN: break
                    sv10, uv10 = calc_sig(nsigv10, nbkgv10, dsigv10, dbkgv10)
                    sv_9_10, _ = sum_z([[sv9,uv9],[sv10,uv10]])
                    if sv_9_10 >= svmax_9_10:
                        jv9stop = 0
                        svmax_9_10 = sv_9_10
                        jv9max___ = jv9
                    else:
                        jv9stop += 1
                    if jv9stop == earlystop: break

                if jv9max___ == 0: continue

                jv11stop = 0
                for jv11 in range(jv10, jv12+1 if nvbf >= 6 else vb+1):
                    nsigv11, dsigv11 = hist_integral(h_sig, jv10, jv11-1)
                    nbkgv11, dbkgv11 = hist_integral(h_bkg, jv10, jv11-1)
                    if nvbf >= 6 and nbkgv11 < minN: continue
                    sv11, uv11 = calc_sig(nsigv11, nbkgv11, dsigv11, dbkgv11)
                    nsigv12, dsigv12 = hist_integral(h_sig, jv11, jv12-1)
                    nbkgv12, dbkgv12 = hist_integral(h_bkg, jv11, jv12-1)
                    if nvbf >= 5 and nbkgv12 < minN: break
                    sv12, uv12 = calc_sig(nsigv12, nbkgv12, dsigv12, dbkgv12)
                    sv_11_12, _ = sum_z([[sv11,uv11],[sv12,uv12]])
                    if sv_11_12 >= svmax_11_12:
                        jv11stop = 0
                        svmax_11_12 = sv_11_12
                        jv11max___ = jv11
                    else:
                        jv11stop += 1
                    if jv11stop == earlystop: break

                if jv11max___ ==0: break

                sv_9_12 = sqrt(svmax_9_10**2 + svmax_11_12**2)
                if sv_9_12 >= svmax_9_12:
                    jv10stop = 0
                    svmax_9_12 = sv_9_12
                    jv9max__, jv10max__, jv11max__ = jv9max___, jv10, jv11max___
                else:
                    jv10stop += 1
                if jv10stop == earlystop: break

            if jv10max__ == 0: continue

            jv14stop = 0
            for jv14 in range(jv12, nscanvbf+1 if nvbf >= 3 else vb+1):

                jv13max___, jv15max___, svmax_13_14, svmax_15_16 = 0, 0, 0, 0
                jv13stop = 0
                for jv13 in range(jv12, jv14+1 if nvbf >= 4 else vb+1):
                    nsigv13, dsigv13 = hist_integral(h_sig, jv12, jv13-1)
                    nbkgv13, dbkgv13 = hist_integral(h_bkg, jv12, jv13-1)
                    if nvbf >= 4 and nbkgv13 < minN: continue
                    sv13, uv13 = calc_sig(nsigv13, nbkgv13, dsigv13, dbkgv13)
                    nsigv14, dsigv14 = hist_integral(h_sig, jv13, jv14-1)
                    nbkgv14, dbkgv14 = hist_integral(h_bkg, jv13, jv14-1)
                    if nvbf >= 3 and nbkgv14 < minN: break
                    sv14, uv14 = calc_sig(nsigv14, nbkgv14, dsigv14, dbkgv14)
                    sv_13_14, _ = sum_z([[sv13,uv13],[sv14,uv14]])
                    if sv_13_14 >= svmax_13_14:
                        jv13stop = 0
                        svmax_13_14 = sv_13_14
                        jv13max___ = jv13
                    else:
                        jv13stop += 1
                    if jv13stop == earlystop: break

                if jv13max___ == 0: continue

                jv15stop = 0
                for jv15 in range(jv14, nscanvbf+1 if nvbf >= 2 else vb+1):
                    nsigv15, dsigv15 = hist_integral(h_sig, jv14, jv15-1)
                    nbkgv15, dbkgv15 = hist_integral(h_bkg, jv14, jv15-1)
                    if nvbf >= 2 and nbkgv15 < minN: continue
                    sv15, uv15 = calc_sig(nsigv15, nbkgv15, dsigv15, dbkgv15)
                    nsigv16, dsigv16 = hist_integral(h_sig, jv15, nscanvbf)
                    nbkgv16, dbkgv16 = hist_integral(h_bkg, jv15, nscanvbf)
                    if nbkgv16 < minN: break
                    sv16, uv16 = calc_sig(nsigv16, nbkgv16, dsigv16, dbkgv16)
                    sv_15_16, _ = sum_z([[sv15,uv15],[sv16,uv16]])
                    if sv_15_16 >= svmax_15_16:
                        jv15stop = 0
                        svmax_15_16 = sv_15_16
                        jv15max___ = jv15
                    else:
                        jv15stop += 1
                    if jv15stop == earlystop: break

                if jv15max___ == 0: break

                sv_13_16 = sqrt(svmax_13_14**2 + svmax_15_16**2)
                if sv_13_16 >= svmax_13_16:
                    jv14stop = 0
                    svmax_13_16 = sv_13_16
                    jv13max__, jv14max__, jv15max__ = jv13max___, jv14, jv15max___
                else:
                    jv14stop+=1
                if jv14stop == earlystop: break

            if jv14max__ == 0: break

            sv_9_16 = sqrt(svmax_9_12**2+svmax_13_16**2)
            if sv_9_16 >= svmax_9_16:
                jv12stop = 0
                svmax_9_16 = sv_9_16
                jv9max_, jv10max_, jv11max_, jv12max_, jv13max_, jv14max_, jv15max_ = jv9max__, jv10max__, jv11max__, jv12, jv13max__, jv14max__, jv15max__
            else:
                jv12stop += 1
            if jv12stop == earlystop: break

        if jv12max_ == 0: break

        sv = sqrt(svmax_1_8**2 + svmax_9_16**2)
        if sv >= svmax:
            jv8stop = 0
            svmax=sv
            jv1max, jv2max, jv3max, jv4max, jv5max, jv6max, jv7max, jv8max, jv9max, jv10max, jv11max, jv12max, jv13max, jv14max, jv15max = jv1max_, jv2max_, jv3max_, jv4max_, jv5max_, jv6max_, jv7max_, jv8, jv9max_, jv10max_, jv11max_, jv12max_, jv13max_, jv14max_, jv15max_
        else:
            jv8stop += 1
        if jv8stop == earlystop: break

    elapsed_time = time.time() - start_time
    tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print '%s |%s| 100/100'%(tt, '>'*100)

    boundaries_VBF = [jv1max, jv2max, jv3max, jv4max, jv5max, jv6max, jv7max, jv8max, jv9max, jv10max, jv11max, jv12max, jv13max, jv14max, jv15max]

    # remove overlapping boundaries
    i=0
    #while boundaries_VBF[0] == vb:
    #    boundaries_VBF.pop(0)
    while i < len(boundaries_VBF) - 1:
        if boundaries_VBF[i+1] == boundaries_VBF[i]:
            boundaries_VBF.pop(i+1)
        else: i += 1

    boundaries_VBF_values = [(i-1.)/nscanvbf for i in boundaries_VBF]

    print '========================================================================='
    print 'Fold number %d' %fold
    print 'The maximal VBF category significance:  %f' %(svmax)
    print 'VBF boundaries: ', boundaries_VBF_values

    #=================================================
    # categorization for general higgs categories

    del h_sig
    del h_bkg
    del h_bkg_left
    del h_bkg_right

    h_sig=TH1F('h_sig','h_sig',nscan,0,1)
    #h_bkg=TH1F('h_bkg','h_bkg',nscan,0,1)

    h_sig.Sumw2()
    #h_bkg.Sumw2()

    t_sig.Draw("bdt_score%s>>h_sig"%('_t' if transform else ''),"weight*%f*((m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,n_fold,fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf))
    #t_bkg.Draw("bdt_score%s>>h_bkg"%('_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.),n_fold,fold, '_t' if transform else '', (vb-1.)/nscanvbf))


    fitboundary = 0.5
    nbin_right = int(nscan*(1-fitboundary))
    nbin_left = int(nscan*fitboundary)
    fitbin = nbin_right  #100

    h_bkg_right = TH1F('h_bkg_right', 'h_bkg_right', fitbin, fitboundary, 1.)
    h_bkg_left = TH1F('h_bkg_left', 'h_bkg_left', nbin_left, 0., fitboundary)
    h_bkg_right.Sumw2()
    h_bkg_left.Sumw2()
    t_bkg.Draw("bdt_score%s>>h_bkg_right"%('_t' if transform else ''),"weight*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f&&bdt_score%s>=%f))"%(n_fold,fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf, '_t' if transform else '', fitboundary))
    t_bkg.Draw("bdt_score%s>>h_bkg_left"%('_t' if transform else ''),"weight*%f*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&(bdt_score_t<%f&&bdt_score_t>=0)&&(eventNumber%%%d!=%d)&&(bdt_score_VBF%s<%f))"%(n_fold/(n_fold-1.) if n_fold != 1 else 1,fitboundary, n_fold,fold if n_fold != 1 else 1, '_t' if transform else '', (vb-1.)/nscanvbf))


    sf = (n_fold/(n_fold-1.) if n_fold != 1 else 1)*0.2723
    #h_bkg_right = fit_BDT(h_bkg_right, fitboundary, fitbin, 'Epoly2', nbin_right, sf)
    h_bkg_right.Scale(sf)

    h_bkg = TH1F('h_bkg','h_bkg',nscan,0.,1.)
    h_bkg.Sumw2()

    h_bkg_list = TList()
    h_bkg_list.Add(h_bkg_left)
    h_bkg_list.Add(h_bkg_right)
    #h_bkg_list.Print()
    h_bkg.Merge(h_bkg_list)




    print 'INFO: scanning all of the possibilities...'
    start_time = time.time()

    j1max, j2max, j3max, j4max, j5max, j6max, j7max, j8max, j9max, j10max, j11max, j12max, j13max, j14max, j15max, smax = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    j8stop = 0
    for j8 in range(1, nscan+1 if nbin >= 9 else 2):

        elapsed_time = time.time() - start_time
        tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print '%s |%s%s| %d/100'%(tt, '>'*(j8*100/nscan), ' '*(100-j8*100/nscan), j8*100/nscan)
        sys.stdout.write("\033[F") # Cursor up one line

        j1max_, j2max_, j3max_, j4max_, j5max_, j6max_, j7max_, smax_1_8 = 0, 0, 0, 0, 0, 0, 0, 0
        j4stop = 0
        for j4 in range(1, j8+1 if nbin >= 13 else 2):

            j1max__, j2max__, j3max__, j5max__, j6max__, j7max__, smax_1_4, smax_5_8 = 0, 0, 0, 0, 0, 0, 0, 0
            j2stop = 0
            for j2 in range(1, j8+1 if nbin >= 15 else 2):

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

    #totalbin=args.nbin+args.vbf
    if not os.path.isdir('cat_opt/%s/%d_%d'%(region, args.vbf, args.nbin)):
        print 'INFO: Creating directory: cat_opt/%s/%d_%d' %(region, args.vbf, args.nbin)
        os.makedirs('cat_opt/%s/%d_%d'%(region, args.vbf, args.nbin))

    with open('cat_opt/%s/%d_%d/%d.json' % (region, args.vbf, args.nbin, args.vb), 'w') as json_file:
        json.dump(outs, json_file)
   

    return



if __name__ == '__main__':
    main()
