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


def apply_weight(categories,region,tsfs,tsfs_VBF=[], VBF=True):

    models=[]
    for i in range(4):
        model = xgb.Booster()
        model.load_model('models/%s_%d.h5'%(region, i))
        models.append(model)
        del model

    if VBF:
        models_VBF=[]
        for i in range(4):
            model = xgb.Booster()
            model.load_model('models/VBF_%s_%d.h5'%(region, i))
            models_VBF.append(model)
            del model

    # loading
    badfile=False
    for category in categories:
        print 'Reading from %s samples...' %(category)
        AppInputsdir="AppInputs/"+category
        npydir='arrays/'+category
        for appinput in os.listdir(AppInputsdir):
            if not appinput.endswith('_%s.root' % (region)):
                continue
            jets=np.load(npydir+'/'+appinput.replace('.root','_-1_jets.npy'))
            if [jets.shape[1],jets.shape[2]]!=[2,5]: 
                print 'ERROR: Dimension of jets npy is not correct!!'
                quit() 
            met=np.load(npydir+'/'+appinput.replace('.root','_-1_met.npy'))
            if [met.shape[1],met.shape[2]]!=[1,2]:
                print 'ERROR: Dimension of met npy is not correct!!'
                quit()
            dijet=np.load(npydir+'/'+appinput.replace('.root','_-1_dijet.npy'))
            if [dijet.shape[1],dijet.shape[2]]!=[1,4]:
                print 'ERROR: Dimension of dijet npy is not correct!!'
                quit()
            dimuon=np.load(npydir+'/'+appinput.replace('.root','_-1_dimuon.npy'))
            if [dimuon.shape[1],dimuon.shape[2]]!=[1,4]:
                print 'ERROR: Dimension of dimuon npy is not correct!!'
                quit()
       
            # get BDT scores for given arrays
            #print 'Getting BDT scores from the trained model...'
            events=[]
            events.append(jets.reshape((jets.shape[0], -1)))
            events.append(met.reshape((met.shape[0], -1)))
            events.append(dijet.reshape((dijet.shape[0], -1)))
            events.append(dimuon.reshape((dimuon.shape[0], -1)))
            events = np.concatenate(events, axis=1)
            dEvents = xgb.DMatrix(events)
 
            scores = []
            scores_t = []
            if VBF:
                scores_VBF = []
                scores_VBF_t = []
            for i in range(4):
                score = models[i].predict(dEvents)
                scores_t.append(tsfs[i].transform(score.reshape(-1,1)).reshape(-1))
                scores.append(score)

                if VBF:
                    score_VBF = models_VBF[i].predict(dEvents)
                    scores_VBF_t.append(tsfs_VBF[i].transform(score_VBF.reshape(-1,1)).reshape(-1))
                    scores_VBF.append(score_VBF)

            infile = ROOT.TFile.Open(AppInputsdir+'/'+appinput)
            intree = infile.Get('AppInput')

            eventNumber = array('L', [0])
            m_mumu = array('f', [0])
            njet = array('I', [0])

            metFinalTrk = array('f', [0])

            Muons_MCPExpReso_Smear_Minv_MuMu_Sigma = array('f', [0])
            Muons_Minv_MuMu_Sigma = array('f', [0])
            Muons_CosThetaStar = array('f', [0])

            Z_PT_FSR = array('f', [0])
            Z_Y_FSR = array('f', [0])
            #Z_Phi_FSR = array('f', [0])

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

            weight=array('d',[0])
            intree.SetBranchAddress('eventNumber', eventNumber)
            intree.SetBranchAddress('m_mumu', m_mumu)
            intree.SetBranchAddress('njet', njet)

            intree.SetBranchAddress('metFinalTrk', metFinalTrk)

            intree.SetBranchAddress('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma)
            intree.SetBranchAddress('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma)
            intree.SetBranchAddress('Muons_CosThetaStar', Muons_CosThetaStar)

            intree.SetBranchAddress('Z_PT_FSR', Z_PT_FSR)
            intree.SetBranchAddress('Z_Y_FSR', Z_Y_FSR)
            #intree.SetBranchAddress('Z_Phi_FSR', Z_Phi_FSR)

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

            # output tree
            bdt_score = array('f', [0])
            if VBF: bdt_score_VBF = array('f', [0])
            bdt_score_t = array('f', [0])
            if VBF: bdt_score_VBF_t = array('f', [0])
            outtree = ROOT.TTree('test', 'test')
            outtree.SetDirectory(0)
            outtree.Branch('eventNumber', eventNumber,'eventNumber/l')
            outtree.Branch('m_mumu', m_mumu,'m_mumu/F')
            outtree.Branch('njet', njet, 'njet/i')

            outtree.Branch('metFinalTrk', metFinalTrk,'metFinalTrk/F')

            outtree.Branch('Z_PT_FSR', Z_PT_FSR, 'Z_PT_FSR/F')
            outtree.Branch('Z_Y_FSR', Z_Y_FSR, 'Z_Y_FSR/F')
            #outtree.Branch('Z_Phi_FSR', Z_Phi, 'Z_Phi/F')

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
            outtree.Branch('bdt_score', bdt_score, 'bdt_score/F')
            if VBF: outtree.Branch('bdt_score_VBF', bdt_score_VBF, 'bdt_score_VBF/F')
            outtree.Branch('bdt_score_t', bdt_score_t, 'bdt_score_t/F')
            if VBF: outtree.Branch('bdt_score_VBF_t', bdt_score_VBF_t, 'bdt_score_VBF_t/F')



            #print 'writing bdt scores into new rootfiles...'
            for i in range(intree.GetEntries()):
                intree.GetEntry(i)
                tag = eventNumber[0]%4
                tag_val = (tag - 1)%4
                bdt_score[0] = 1.0 - scores[tag][i]
                bdt_score_t[0] = 1 - scores_t[tag][i]
                if VBF:
                    bdt_score_VBF[0] = 1.0 - scores_VBF[tag][i]
                    bdt_score_VBF_t[0] = 1 - scores_VBF_t[tag][i]
                outtree.Fill()

            if not os.path.isdir('outputs'):
                print 'INFO: Creating new folder: \"outputs\"'
                os.makedirs("outputs")
            if not os.path.isdir('outputs/model_%s' %(region)):
                print 'INFO: Creating model container:  model_%s' %(region)
                os.makedirs('outputs/model_%s' %(region))
            if not os.path.isdir('outputs/model_%s/%s' %(region,category)):
                print 'INFO: Creating in the model container a new category:  %s' %(category)
                os.makedirs('outputs/model_%s/%s' %(region,category))


            outfile = ROOT.TFile('outputs/model_%s/%s/%s' % (region,category,appinput), 'RECREATE')
            outtree.Write()
            outfile.Close()

            infile.Close()

        print 'Combining all the root-files in the category:  ', category
        if os.path.isdir('outputs/model_%s/%s'%(region,category)):
            os.system("hadd -f outputs/model_%s/%s.root outputs/model_%s/%s/*.root" %(region,category,region,category))

    return



def main():
    args=getArgs()
    

    sigs = ['ggF','VBF','VH','ttH']
    bkgs = ['data']

    region = args.region

    if region in ['all_jet', 'two_jet']:
        VBF = True
    else:
        VBF = False

    print '=============================================================================='
    print 'INFO:  Using the model: %s.h'%(region)
    if VBF: 
        print 'INFO:  Using the model: VBF_%s.h'%(region) 
    print '------------------------------------------------------------------------------'

    if VBF:
        tsfs_VBF = []
        for i in range(4):
            tsf = pickle.load(open('models/score_transformer_%s_VBF_%d.pkl'%(region, i), "rb" ) )
            tsfs_VBF.append(tsf)

    tsfs = []
    for i in range(4):
        tsf = pickle.load(open('models/score_transformer_%s_%d.pkl'%(region, i), "rb" ) )
        tsfs.append(tsf)

    if VBF:
        apply_weight(sigs+bkgs, region, tsfs, tsfs_VBF=tsfs_VBF)
    else:
        apply_weight(sigs+bkgs, region, tsfs, VBF=False)

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
