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
from math import sqrt, log
import math
import sys
import time
import json

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'two_jet', help = 'Region to process')
    parser.add_argument('-b', '--nbin', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 3, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
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

def hist_integral(hist,i,j,k,l):
    err = ROOT.Double()
    if i>j or k>l:
        n=0
        err=0
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

def gettingsig(region, boundaries_VBF, boundaries, nscanvbf, nscan, nfold, transform):

    f_bkg = ROOT.TFile('outputs/%s/bkg.root' % (region))
    t_bkg = f_bkg.Get('test')

    h_bkg = []
    for fold in range(nfold):
        h_bkg.append(ROOT.TH2F('h_bkg_%d'%(fold),'h_bkg_%d'%(fold),nscanvbf,0.,1.,nscan,0.,1.))
        h_bkg[fold].Sumw2()
        t_bkg.Draw("bdt_score%s:bdt_score_VBF%s>>h_bkg_%d"%('_t' if transform else '','_t' if transform else '', fold),"weight*1*(0.2723)*((m_mumu>=110&&m_mumu<=180)&&!(m_mumu>=120&&m_mumu<=130)&&eventNumber%%%d==%d)"%(nfold,fold))

    nbkg = []
    for i in range(len(boundaries_VBF[0])):
        nbkg.append(sum_hist_integral([hist_integral(h_bkg[fold], boundaries_VBF[fold][i], boundaries_VBF[fold][i+1]-1 if i != len(boundaries_VBF[0])-1 else nscanvbf+1, 1, nscan+1) for fold in range(nfold)]))
    for i in range(len(boundaries[0])):
        nbkg.append(sum_hist_integral([hist_integral(h_bkg[fold], 1, boundaries_VBF[fold][0]-1, boundaries[fold][i], boundaries[fold][i+1]-1 if i != len(boundaries[0])-1 else nscan+1) for fold in range(nfold)]))


    f_sig = ROOT.TFile('outputs/%s/sig.root' % (region))
    t_sig = f_sig.Get('test')

    h_sig = []
    for fold in range(nfold):
        h_sig.append(ROOT.TH2F('h_sig_%d'%(fold),'h_bkg_%d'%(fold),nscanvbf,0.,1.,nscan,0.,1.))
        h_sig[fold].Sumw2()
        t_sig.Draw("bdt_score%s:bdt_score_VBF%s>>h_sig_%d"%('_t' if transform else '', '_t' if transform else '', fold),"weight*1*(m_mumu>=120&&m_mumu<=130&&eventNumber%%%d==%d)"%(nfold,fold))

    nsig = []
    for i in range(len(boundaries_VBF[0])):
        nsig.append(sum_hist_integral([hist_integral(h_sig[fold], boundaries_VBF[fold][i], boundaries_VBF[fold][i+1]-1 if i != len(boundaries_VBF[0])-1 else nscanvbf+1, 1, nscan+1) for fold in range(nfold)]))
    for i in range(len(boundaries[0])):
        nsig.append(sum_hist_integral([hist_integral(h_sig[fold], 1, boundaries_VBF[fold][0]-1, boundaries[fold][i], boundaries[fold][i+1]-1 if i != len(boundaries[0])-1 else nscan+1) for fold in range(nfold)]))       
      
    zs = [ calc_sig(nsig[i][0], nbkg[i][0], nsig[i][1], nbkg[i][1]) for i in range(len(nsig))]
    s, u = sum_z(zs)
 
    print 'Signal yield: ', nsig
    print 'BKG yield: ', nbkg
    print 'Significance:  %f +- %f'%(s,abs(u))

    return s, abs(u)


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
    
    s, u = gettingsig(region, boundaries_VBF, boundaries, nscanvbf, nscan, nfold, transform)

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
