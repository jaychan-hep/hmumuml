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
    parser.add_argument('-v', '--var', action = 'store', choices = ['Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', 'Muons_Minv_MuMu_Sigma'], default = 'Muons_Minv_MuMu_Sigma', help = 'reso var')
    return parser.parse_args()


def apply_weight(categories, region, boundaries, var):

    nbin = len(boundaries[0])
    print 'Number of bins: ', nbin

    # loading
    for category in categories:

        print 'Reading from %s samples...' %(category)

        infile = ROOT.TFile.Open('outputs_2D/model_%s/%s.root'%(region, category))
        intree = infile.Get('test')
        
        eventNumber = array('L', [0])
        m_mumu = array('f', [0])
        njet = array('I', [0])

        metFinalTrk = array('f', [0])

        Muons_MCPExpReso_Smear_Minv_MuMu_Sigma = array('f', [0])
        Muons_Minv_MuMu_Sigma = array('f', [0])
        Muons_CosThetaStar = array('f', [0])

        Z_PT = array('f', [0])
        Z_Eta = array('f', [0])
        #Z_Phi = array('f', [0])

        Jets_PT_Lead = array('f', [0])
        Jets_PT_Sub = array('f', [0])
        Jets_Eta_Lead = array('f', [0])
        Jets_Eta_Sub = array('f', [0])
        DeltaPhi_mumuj1 = array('f', [0])
        DeltaPhi_mumuj2 = array('f', [0])
        Jets_PT_jj = array('f', [0])
        Jets_Eta_jj = array('f', [0])
        DeltaPhi_mumujj = array('f', [0])
        Jets_Minv_jj = array('f', [0])

        weight = array('d',[0])
        bdt_score_t = array('f', [0])
        if region == 'two_jet': bdt_score_VBF_t = array('f', [0])
        bdt_category = array('f', [0])
        reso_category = array('f', [0])

        intree.SetBranchAddress('eventNumber', eventNumber)
        intree.SetBranchAddress('m_mumu', m_mumu)
        intree.SetBranchAddress('njet', njet)

        intree.SetBranchAddress('metFinalTrk', metFinalTrk)

        intree.SetBranchAddress('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma)
        intree.SetBranchAddress('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma)
        intree.SetBranchAddress('Muons_CosThetaStar', Muons_CosThetaStar)

        intree.SetBranchAddress('Z_PT', Z_PT)
        intree.SetBranchAddress('Z_Eta', Z_Eta)
        #intree.SetBranchAddress('Z_Phi', Z_Phi)

        intree.SetBranchAddress('Jets_PT_Lead', Jets_PT_Lead)
        intree.SetBranchAddress('Jets_PT_Sub', Jets_PT_Sub)
        intree.SetBranchAddress('Jets_Eta_Lead', Jets_Eta_Lead)
        intree.SetBranchAddress('Jets_Eta_Sub', Jets_Eta_Sub)
        intree.SetBranchAddress('DeltaPhi_mumuj1', DeltaPhi_mumuj1)
        intree.SetBranchAddress('DeltaPhi_mumuj2', DeltaPhi_mumuj2)
        intree.SetBranchAddress('Jets_PT_jj', Jets_PT_jj)
        intree.SetBranchAddress('Jets_Eta_jj', Jets_Eta_jj)
        intree.SetBranchAddress('DeltaPhi_mumujj', DeltaPhi_mumujj)
        intree.SetBranchAddress('Jets_Minv_jj', Jets_Minv_jj)

        intree.SetBranchAddress('weight', weight)
        intree.SetBranchAddress('bdt_score_t', bdt_score_t)
        if region == 'two_jet': intree.SetBranchAddress('bdt_score_VBF_t', bdt_score_VBF_t)
        intree.SetBranchAddress('bdt_category', bdt_category)

        # output tree
        outtree = ROOT.TTree('test', 'test')
        outtree.SetDirectory(0)
        outtree.Branch('eventNumber', eventNumber,'eventNumber/l')
        outtree.Branch('m_mumu', m_mumu,'m_mumu/F')
        outtree.Branch('njet', njet, 'njet/i')

        outtree.Branch('metFinalTrk', metFinalTrk,'metFinalTrk/F')

        outtree.Branch('Z_PT', Z_PT, 'Z_PT/F')
        outtree.Branch('Z_Eta', Z_Eta, 'Z_Eta/F')
        #outtree.Branch('Z_Phi', Z_Phi, 'Z_Phi/F')

        outtree.Branch('Jets_PT_Lead', Jets_PT_Lead, 'Jets_PT_Lead/F')
        outtree.Branch('Jets_PT_Sub', Jets_PT_Sub, 'Jets_PT_Sub/F')
        outtree.Branch('Jets_Eta_Lead', Jets_Eta_Lead, 'Jets_Eta_Lead/F')
        outtree.Branch('Jets_Eta_Sub', Jets_Eta_Sub, 'Jets_Eta_Sub/F')
        outtree.Branch('DeltaPhi_mumuj1', DeltaPhi_mumuj1, 'DeltaPhi_mumuj1/F')
        outtree.Branch('DeltaPhi_mumuj2', DeltaPhi_mumuj2, 'DeltaPhi_mumuj2/F')

        outtree.Branch('Jets_PT_jj', Jets_PT_jj, 'Jets_PT_jj/F')
        outtree.Branch('Jets_Eta_jj', Jets_Eta_jj, 'Jets_Eta_jj/F')
        outtree.Branch('DeltaPhi_mumujj', DeltaPhi_mumujj, 'DeltaPhi_mumujj/F')
        outtree.Branch('Jets_Minv_jj', Jets_Minv_jj, 'Jets_Minv_jj/F')

        outtree.Branch('Muons_CosThetaStar', Muons_CosThetaStar, 'Muons_CosThetaStar/F')

        outtree.Branch('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma, 'Muons_MCPExpReso_Smear_Minv_MuMu_Sigma/F')
        outtree.Branch('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma, 'Muons_Minv_MuMu_Sigma/F')

        outtree.Branch('weight', weight, 'weight/D')
        outtree.Branch('bdt_score_t', bdt_score_t, 'bdt_score_t/F')
        if region == 'two_jet': outtree.Branch('bdt_score_VBF_t', bdt_score_VBF_t, 'bdt_score_VBF_t/F')
        outtree.Branch('bdt_category', bdt_category, 'bdt_category/F')
        outtree.Branch('reso_category', reso_category, 'reso_category/F')

        #print 'writing bdt scores into new rootfiles...'
        for i in range(intree.GetEntries()):

            intree.GetEntry(i)
            tag = eventNumber[0]%4

            if var == 'Muons_Minv_MuMu_Sigma': value = Muons_Minv_MuMu_Sigma[0]/m_mumu[0]
            elif var == 'Muons_MCPExpReso_Smear_Minv_MuMu_Sigma': value = Muons_MCPExpReso_Smear_Minv_MuMu_Sigma[0]/m_mumu[0]

            cat = int(bdt_category[0])
            if boundaries[0][cat-1] == 0: reso_category[0] = 1
            elif value >= boundaries[tag][cat-1]: reso_category[0] = 0
            else: reso_category[0] = 1
         
            outtree.Fill()

        if not os.path.isdir('outputs_reso'):
            print 'INFO: Creating new folder: \"outputs_reso\"'
            os.makedirs("outputs_reso")
        if not os.path.isdir('outputs_reso/model_%s' %(region)):
            print 'INFO: Creating folder:  model_%s' %(region)
            os.makedirs('outputs_reso/model_%s' %(region))

        outfile = ROOT.TFile('outputs_reso/model_%s/%s.root' % (region, category), 'RECREATE')
        outtree.Write()
        outfile.Close()

        infile.Close()

    return



def main():
    args=getArgs()
    

    sigs = ['ggF','VBF','VH','ttH']
    bkgs = ['data']

    region = args.region

    if region == 'zero_jet': boundaries = [[0.01722, 0.01682, 0.018000000000000002, 0.01846, 0.01998], [0.01726, 0.01668, 0.01866, 0.01864, 0.02064], [0.01722, 0.01726, 0.018000000000000002, 0.01944, 0.020540000000000003], [0.01726, 0.017320000000000002, 0.01794, 0.01864, 0.019700000000000002]]
    elif region == 'one_jet': boundaries = [[0.01718, 0.01942, 0.01842, 0.01784, 0.02048], [0.021859999999999997, 0.01816, 0.01996, 0.01806, 0.018840000000000003], [0.01718, 0.018840000000000003, 0.01998, 0.0176, 0.020540000000000003], [0.01718, 0.018840000000000003, 0.01996, 0.018500000000000003, 0.020540000000000003]]
    elif region == 'two_jet': boundaries = [[0, 0.02038, 0.0187, 0.01944, 0.01752, 0.02098, 0.02098, 0.02096, 0.01962, 0.02064], [0, 0.02058, 0.0187, 0.01944, 0.01764, 0.02052, 0.02108, 0.02038, 0.01884, 0.02106], [0, 0.02038, 0.01874, 0.01948, 0.01752, 0.02312, 0.02256, 0.0209, 0.01962, 0.0216], [0, 0.02038, 0.01898, 0.01952, 0.01752, 0.02322, 0.02054, 0.02092, 0.01884, 0.0216]]


    apply_weight(sigs+bkgs, region, boundaries, args.var)

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
