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
import json

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['zero_jet','one_jet','two_jet','all_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-b', '--nbin', type = int, default = 4, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 0, choices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    return parser.parse_args()


def apply_weight(categories, region, boundaries, boundaries_v):

    nbin = len(boundaries[0]) + len(boundaries_v[0]) 
    nVBF = len(boundaries_v[0])
    nggF = len(boundaries[0]) 
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
        EventWeight_MCEventWeight = array('f', [0])
        bdt_score_t = array('f', [0])
        if nVBF != 0: bdt_score_VBF_t = array('f', [0])
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
        intree.SetBranchAddress('EventWeight_MCEventWeight', EventWeight_MCEventWeight)
        intree.SetBranchAddress('bdt_score_t', bdt_score_t)
        if nVBF != 0: intree.SetBranchAddress('bdt_score_VBF_t', bdt_score_VBF_t)

        # outfile
        if not os.path.isdir('outputs_2D'):
            print 'INFO: Creating new folder: \"outputs_2D\"'
            os.makedirs("outputs_2D")
        if not os.path.isdir('outputs_2D/model_%s' %(region)):
            print 'INFO: Creating folder:  model_%s' %(region)
            os.makedirs('outputs_2D/model_%s' %(region))

        outfile = ROOT.TFile('outputs_2D/model_%s/%s.root' % (region, category), 'RECREATE')

        # output tree
        outtree = ROOT.TTree('test', 'test')
        outtree.SetAutoSave(10000)
        #outtree.SetDirectory(0)
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
        outtree.Branch('EventWeight_MCEventWeight', EventWeight_MCEventWeight, 'EventWeight_MCEventWeight/F')
        outtree.Branch('bdt_score_t', bdt_score_t, 'bdt_score_t/F')
        if nVBF != 0: outtree.Branch('bdt_score_VBF_t', bdt_score_VBF_t, 'bdt_score_VBF_t/F')
        outtree.Branch('bdt_category', bdt_category, 'bdt_category/F')

        #print 'writing bdt scores into new rootfiles...'
        for i in range(intree.GetEntries()):

            intree.GetEntry(i)
            tag = eventNumber[0]%len(boundaries)

            if nVBF == 0 or bdt_score_VBF_t[0] < float(boundaries_v[tag][0]):

                if bdt_score_t[0] < boundaries[tag][0]: continue
                bdt_category[0] = nVBF + 1
                for j in range(1, nggF):
                    if bdt_score_t[0] < float(boundaries[tag][j]):
                        bdt_category[0] = nbin - j + 1
                        break
            else:
                bdt_category[0] = 1
                for j in range(nVBF-1):
                    if bdt_score_VBF_t[0] < float(boundaries_v[tag][j+1]):
                        bdt_category[0] = nVBF - j
                        break
         
            outtree.Fill()

        outtree.Write()
        outfile.Close()

        infile.Close()

    os.system('cd outputs_2D/model_%s; hadd sig.root ggF.root VBF.root VH.root ttH.root'%(region))

    return



def main():
    args=getArgs()
    

    sigs = ['ggF','VBF','VH','ttH']
    bkgs = ['data']
    bkgs += ['Z', 'ttbar', 'diboson', 'stop']

    region = args.region

    nvbf, nggf = args.vbf, args.nbin
    with open('significances/%s/%s%d.json'%(region, '%d_'%nvbf if nvbf != 0 else '', nggf)) as f:
        data = json.load(f)
    #transform = data['transform']
    #nfold = data['nfold']
    #nscanvbf = data['nscanvbf']
    #nscan = data['nscan']
    #floatB = data['floatB']
    #minN = data['minN']
    boundaries_v = data['boundaries_VBF_values'] if nvbf != 0 else [[]]
    boundaries = data['boundaries_values']
    
    apply_weight(sigs+bkgs, region, boundaries, boundaries_v)

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
