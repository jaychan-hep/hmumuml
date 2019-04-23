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
    parser.add_argument('-r', '--region', action='store', choices=['zero_jet','one_jet','two_jet','all_jet'], default='all_jet', help='Region to process')
    return parser.parse_args()


def apply_weight(categories, region, boundaries):

    nbin = len(boundaries[0])
    print 'Number of bins: ', nbin

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
        #intree.SetBranchAddress('bdt_score_VBF_t', bdt_score_VBF_t)

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
        #outtree.Branch('bdt_score_VBF_t', bdt_score_VBF_t, 'bdt_score_VBF_t/F')
        outtree.Branch('bdt_category', bdt_category, 'bdt_category/F')

        #print 'writing bdt scores into new rootfiles...'
        for i in range(intree.GetEntries()):

            intree.GetEntry(i)
            tag = eventNumber[0]%4

            if bdt_score_t[0] < boundaries[tag][0]: continue
            bdt_category[0] = 1
            for j in range(1, nbin):
                if bdt_score_t[0] < float(boundaries[tag][j]):
                    bdt_category[0] = nbin - j + 1
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

    if region == 'zero_jet':
        # 6 bins
        #b0=[0.0, 0.14, 0.36, 0.58, 0.77, 0.91]
        #b1=[0.0, 0.13, 0.36, 0.61, 0.78, 0.92]
        #b2=[0.0, 0.13, 0.33, 0.54, 0.75, 0.92]
        #b3=[0.0, 0.13, 0.34, 0.54, 0.75, 0.92]

        # 5 bins
        #b0=[0.0, 0.16, 0.4, 0.66, 0.86]
        #b1=[0.0, 0.14, 0.38, 0.65, 0.87]
        #b2=[0.0, 0.16, 0.4, 0.66, 0.87]
        #b3=[0.0, 0.15, 0.4, 0.65, 0.86]

        # 4 bins
        #b0=[0.0, 0.24, 0.56, 0.85]
        #b1=[0.0, 0.24, 0.56, 0.85]
        #b2=[0.0, 0.23, 0.54, 0.85]
        #b3=[0.0, 0.24, 0.54, 0.85]

        # 3 bins
        #b0=[0.0, 0.36, 0.77]
        #b1=[0.0, 0.36, 0.75]
        #b2=[0.0, 0.34, 0.75]
        #b3=[0.0, 0.34, 0.75]

        # 2 bins
        b0=[0.0, 0.52]
        b1=[0.0, 0.52]
        b2=[0.0, 0.54]
        b3=[0.0, 0.53]

    elif region == 'one_jet':
        # 6 bins
        #b0=[0.0, 0.12, 0.38, 0.6, 0.77, 0.91]
        #b1=[0.0, 0.12, 0.38, 0.61, 0.78, 0.92]
        #b2=[0.0, 0.13, 0.38, 0.61, 0.78, 0.93]
        #b3=[0.0, 0.12, 0.38, 0.61, 0.79, 0.93]

        # 5 bins
        #b0=[0.0, 0.2, 0.53, 0.76, 0.91]
        #b1=[0.0, 0.2, 0.54, 0.76, 0.92]
        #b2=[0.0, 0.2, 0.52, 0.76, 0.92]
        #b3=[0.0, 0.2, 0.54, 0.77, 0.93]

        # 4 bins
        #b0=[0.0, 0.27, 0.62, 0.86]
        #b1=[0.0, 0.26, 0.6, 0.84]
        #b2=[0.0, 0.27, 0.61, 0.87]
        #b3=[0.0, 0.27, 0.61, 0.87]

        # 3 bins
        #b0=[0.0, 0.38, 0.77]
        #b1=[0.0, 0.38, 0.78]
        #b2=[0.0, 0.38, 0.78]
        #b3=[0.0, 0.38, 0.78]

        # 2 bins
        b0=[0.0, 0.62]
        b1=[0.0, 0.6]
        b2=[0.0, 0.61]
        b3=[0.0, 0.61]

    elif region == 'two_jet':
        b0=[]
        b1=[]
        b2=[]
        b3=[]

    boundaries=[]
    boundaries.append(b0)
    boundaries.append(b1)
    boundaries.append(b2)
    boundaries.append(b3)

    apply_weight(sigs+bkgs, region, boundaries)

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
