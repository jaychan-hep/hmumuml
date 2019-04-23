#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.12.2018
#
#  Modified from ttHyy XGBosst anaylsis codes
#
#
#
#
import os
from argparse import ArgumentParser
import math
import ROOT
from root_numpy import root2array
from array import array
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from keras import models
import xgboost as xgb
import matplotlib.pyplot as plt


def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['zero_jet','one_jet','two_jet','all_jet'], default='two_jet', help='Region to process')
    return parser.parse_args()


def apply_weight(categories, region, boundaries, boundaries_v):

    nbin = len(boundaries[0]) + len(boundaries_v[0]) + 1
    nVBF = len(boundaries_v[0])
    nggF = len(boundaries[0]) + 1
    print 'Number of VBF bins: ', nVBF
    print 'Number of ggF bins: ', nggF

    # loading
    for category in categories:

        print 'Reading from %s samples...' %(category)

        infile = ROOT.TFile.Open('outputs/model_%s/%s.root'%(region, category))
        intree = infile.Get('test')
        
        eventNumber = array('L', [0])
        m_mumu = array('f', [0])
        njet = array('I', [0])

        metFinalTrk = array('f', [0])

        Muons_MCPExpReso_Smear_Minv_MuMu_Sigma = array('f', [0])
        Muons_Minv_MuMu_Sigma = array('f', [0])
        Muons_CosThetaStar = array('f', [0])

        Z_PT_FSR = array('f', [0])
        Z_Y_FSR = array('f', [0])
        #Z_Phi = array('f', [0])

        Jets_PT_Lead = array('f', [0])
        Jets_PT_Sub = array('f', [0])
        Jets_Eta_Lead = array('f', [0])
        Jets_Eta_Sub = array('f', [0])
        DeltaPhi_mumuj1 = array('f', [0])
        DeltaPhi_mumuj2 = array('f', [0])
        Jets_PT_jj = array('f', [0])
        Jets_Y_jj = array('f', [0])
        DeltaPhi_mumujj = array('f', [0])
        Jets_Minv_jj = array('f', [0])

        weight = array('d',[0])
        bdt_score_t = array('f', [0])
        bdt_score_VBF_t = array('f', [0])
        bdt_category = array('f', [0])

        intree.SetBranchAddress('eventNumber', eventNumber)
        intree.SetBranchAddress('m_mumu', m_mumu)
        intree.SetBranchAddress('njet', njet)

        intree.SetBranchAddress('metFinalTrk', metFinalTrk)

        intree.SetBranchAddress('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma)
        intree.SetBranchAddress('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma)
        intree.SetBranchAddress('Muons_CosThetaStar', Muons_CosThetaStar)

        intree.SetBranchAddress('Z_PT_FSR', Z_PT_FSR)
        intree.SetBranchAddress('Z_Y_FSR', Z_Y_FSR)
        #intree.SetBranchAddress('Z_Phi', Z_Phi)

        intree.SetBranchAddress('Jets_PT_Lead', Jets_PT_Lead)
        intree.SetBranchAddress('Jets_PT_Sub', Jets_PT_Sub)
        intree.SetBranchAddress('Jets_Eta_Lead', Jets_Eta_Lead)
        intree.SetBranchAddress('Jets_Eta_Sub', Jets_Eta_Sub)
        intree.SetBranchAddress('DeltaPhi_mumuj1', DeltaPhi_mumuj1)
        intree.SetBranchAddress('DeltaPhi_mumuj2', DeltaPhi_mumuj2)
        intree.SetBranchAddress('Jets_PT_jj', Jets_PT_jj)
        intree.SetBranchAddress('Jets_Y_jj', Jets_Y_jj)
        intree.SetBranchAddress('DeltaPhi_mumujj', DeltaPhi_mumujj)
        intree.SetBranchAddress('Jets_Minv_jj', Jets_Minv_jj)

        intree.SetBranchAddress('weight', weight)
        intree.SetBranchAddress('bdt_score_t', bdt_score_t)
        intree.SetBranchAddress('bdt_score_VBF_t', bdt_score_VBF_t)

        # output tree
        outtree = ROOT.TTree('test', 'test')
        outtree.SetDirectory(0)
        outtree.Branch('eventNumber', eventNumber,'eventNumber/l')
        outtree.Branch('m_mumu', m_mumu,'m_mumu/F')
        outtree.Branch('njet', njet, 'njet/i')

        outtree.Branch('metFinalTrk', metFinalTrk,'metFinalTrk/F')

        outtree.Branch('Z_PT_FSR', Z_PT_FSR, 'Z_PT_FSR/F')
        outtree.Branch('Z_Y_FSR', Z_Y_FSR, 'Z_Y_FSR/F')
        #outtree.Branch('Z_Phi', Z_Phi, 'Z_Phi/F')

        outtree.Branch('Jets_PT_Lead', Jets_PT_Lead, 'Jets_PT_Lead/F')
        outtree.Branch('Jets_PT_Sub', Jets_PT_Sub, 'Jets_PT_Sub/F')
        outtree.Branch('Jets_Eta_Lead', Jets_Eta_Lead, 'Jets_Eta_Lead/F')
        outtree.Branch('Jets_Eta_Sub', Jets_Eta_Sub, 'Jets_Eta_Sub/F')
        outtree.Branch('DeltaPhi_mumuj1', DeltaPhi_mumuj1, 'DeltaPhi_mumuj1/F')
        outtree.Branch('DeltaPhi_mumuj2', DeltaPhi_mumuj2, 'DeltaPhi_mumuj2/F')

        outtree.Branch('Jets_PT_jj', Jets_PT_jj, 'Jets_PT_jj/F')
        outtree.Branch('Jets_Y_jj', Jets_Y_jj, 'Jets_Y_jj/F')
        outtree.Branch('DeltaPhi_mumujj', DeltaPhi_mumujj, 'DeltaPhi_mumujj/F')
        outtree.Branch('Jets_Minv_jj', Jets_Minv_jj, 'Jets_Minv_jj/F')

        outtree.Branch('Muons_CosThetaStar', Muons_CosThetaStar, 'Muons_CosThetaStar/F')

        outtree.Branch('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma, 'Muons_MCPExpReso_Smear_Minv_MuMu_Sigma/F')
        outtree.Branch('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma, 'Muons_Minv_MuMu_Sigma/F')

        outtree.Branch('weight', weight, 'weight/D')
        outtree.Branch('bdt_score_t', bdt_score_t, 'bdt_score_t/F')
        outtree.Branch('bdt_score_VBF_t', bdt_score_VBF_t, 'bdt_score_VBF_t/F')
        outtree.Branch('bdt_category', bdt_category, 'bdt_category/F')

        #print 'writing bdt scores into new rootfiles...'
        for i in range(intree.GetEntries()):

            intree.GetEntry(i)
            tag = eventNumber[0]%4

            if bdt_score_VBF_t[0] < float(boundaries_v[tag][0]):
                bdt_category[0] = nVBF + 1
                for j in range(nggF-1):
                    if bdt_score_t[0] < float(boundaries[tag][j]):
                        bdt_category[0] = nbin - j
                        break
            else:
                bdt_category[0] = 1
                for j in range(nVBF-1):
                    if bdt_score_VBF_t[0] < float(boundaries_v[tag][j+1]):
                        bdt_category[0] = nVBF - j
                        break
         
            outtree.Fill()

        if not os.path.isdir('outputs_2D'):
            print 'INFO: Creating new folder: \"outputs_2D\"'
            os.makedirs("outputs_2D")
        if not os.path.isdir('outputs_2D/model_%s' %(region)):
            print 'INFO: Creating folder:  model_%s' %(region)
            os.makedirs('outputs_2D/model_%s' %(region))

        outfile = ROOT.TFile('outputs_2D/model_%s/%s.root' % (region, category), 'RECREATE')
        outtree.Write()
        outfile.Close()

        infile.Close()

    return



def main():
    args=getArgs()
    

    sigs = ['ggF','VBF','VH','ttH']
    bkgs = ['data']

    region = args.region

    # 6+6
    #b0='0.050 0.210 0.400 0.560 0.700'
    #b1='0.050 0.220 0.400 0.550 0.690'
    #b2='0.060 0.220 0.380 0.530 0.640'
    #b3='0.050 0.200 0.380 0.530 0.660'
    #c0='0.630 0.710 0.770 0.870 0.920 0.950'
    #c1='0.620 0.770 0.880 0.920 0.970 0.980'
    #c2='0.610 0.700 0.770 0.870 0.920 0.950'
    #c3='0.620 0.770 0.830 0.880 0.900 0.920'

    # 5+5
    #b0='0.050 0.210 0.400 0.570'
    #b1='0.050 0.220 0.400 0.560'
    #b2='0.060 0.220 0.400 0.570'
    #b3='0.080 0.280 0.470 0.650'
    #c0='0.630 0.770 0.870 0.920 0.950'
    #c1='0.610 0.770 0.880 0.920 0.970'
    #c2='0.620 0.770 0.870 0.920 0.950'
    #c3='0.620 0.770 0.880 0.900 0.920'

    # 4+4
    #b0='0.110 0.340 0.560'
    #b1='0.110 0.340 0.560'
    #b2='0.110 0.340 0.560'
    #b3='0.110 0.320 0.560'
    #c0='0.630 0.770 0.870 0.920'
    #c1='0.610 0.770 0.920 0.970'
    #c2='0.620 0.770 0.870 0.920'
    #c3='0.620 0.770 0.870 0.920'

    # 3+3
    #b0='0.200 0.480'
    #b1='0.190 0.470'
    #b2='0.200 0.480'
    #b3='0.200 0.470'
    #c0='0.620 0.770 0.920'
    #c1='0.600 0.770 0.920'
    #c2='0.610 0.770 0.920'
    #c3='0.610 0.770 0.920'

    # 2+2
    b0='0.34'
    b1='0.34'
    b2='0.38'
    b3='0.32'
    c0='0.670 0.920'
    c1='0.620 0.880'
    c2='0.700 0.920'
    c3='0.620 0.880'

    # minN=10

    boundaries=[]
    boundaries.append(b0.split())
    boundaries.append(b1.split())
    boundaries.append(b2.split())
    boundaries.append(b3.split())

    boundaries_v=[]
    boundaries_v.append(c0.split())
    boundaries_v.append(c1.split())
    boundaries_v.append(c2.split())
    boundaries_v.append(c3.split())

    apply_weight(sigs+bkgs, region, boundaries, boundaries_v)

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
