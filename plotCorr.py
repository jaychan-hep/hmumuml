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
#from root_numpy import root2array
from array import array
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
#from keras import models
import xgboost as xgb

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['zero_jet','one_jet','two_jet','all_jet'], default='all_jet', help='Region to process')
    return parser.parse_args()


def apply_weight(categories, region, title):

    data = []

    # loading
    for category in categories:

        print 'Reading from %s samples...' %(category)

        outputsdir="outputs/model_%s/"%(region)+category
        npydir='arrays/'+category
        for output in os.listdir(outputsdir):
            if not output.endswith('_%s.root' % (region)):
                continue
            jets=np.load(npydir+'/'+output.replace('.root','_-1_jets.npy'))
            if [jets.shape[1],jets.shape[2]]!=[2,3]:
                print 'ERROR: Dimension of jets npy is not correct!!'
                error
            met=np.load(npydir+'/'+output.replace('.root','_-1_met.npy'))
            if [met.shape[1],met.shape[2]]!=[1,2]:
                print 'ERROR: Dimension of met npy is not correct!!'
                error
            #muons=np.load(npydir+'/'+appinput.replace('.root','_-1_muons.npy'))
            #if [muons.shape[1],muons.shape[2]]!=[2,4]:
            #    print 'ERROR: Dimension of muons npy is not correct!!'
            #    error
            dijet=np.load(npydir+'/'+output.replace('.root','_-1_dijet.npy'))
            if [dijet.shape[1],dijet.shape[2]]!=[1,4]:
                print 'ERROR: Dimension of dijet npy is not correct!!'
                error
            dimuon=np.load(npydir+'/'+output.replace('.root','_-1_dimuon.npy'))
            if [dimuon.shape[1],dimuon.shape[2]]!=[1,4]:
                print 'ERROR: Dimension of dimuon npy is not correct!!'
                error
            #extras=np.load(npydir+'/'+output.replace('.root','_-1_extras.npy'))
            #if [extras.shape[1],extras.shape[2]]!=[1,6]:
            #    print 'ERROR: Dimension of extras npy is not correct!!'
            #    error

            events=[]
            events.append(jets.reshape((jets.shape[0], -1)))
            events.append(met.reshape((met.shape[0], -1)))
            #events.append(muons.reshape((muons.shape[0], -1)))
            events.append(dijet.reshape((dijet.shape[0], -1)))
            events.append(dimuon.reshape((dimuon.shape[0], -1)))
            #events.append(extras.reshape((extras.shape[0], -1)))
            events = np.concatenate(events, axis=1)


            infile = ROOT.TFile.Open('outputs/model_%s/%s/%s'%(region, category, output))
            intree = infile.Get('test')
        
            m_mumu = array('f', [0])
            bdt_score_t = array('f', [0])
            bdt_score_VBF_t = array('f', [0])

            intree.SetBranchAddress('m_mumu', m_mumu)
            intree.SetBranchAddress('bdt_score_t', bdt_score_t)
            if region == 'two_jet': intree.SetBranchAddress('bdt_score_VBF_t', bdt_score_VBF_t)

            scores = []
            scores_VBF = []
            mass = []

            #print 'writing bdt scores into new rootfiles...'
            for i in range(intree.GetEntries()):

                intree.GetEntry(i)

                scores.append([bdt_score_t[0]])
                if region == 'two_jet': scores_VBF.append([bdt_score_VBF_t[0]])
                mass.append([m_mumu[0]])

            infile.Close()

            if region == 'two_jet': events = np.concatenate([events, scores, scores_VBF, mass], axis=1)
            else: events = np.concatenate([events, scores, mass], axis=1)
            data.append(events)

    data = np.concatenate(data)

    if region == 'two_jet': dummy_list = [7, 14]
    elif region == 'one_jet': dummy_list = [3, 4, 5, 6, 7, 8, 9, 10, 11, 14]
    elif region == 'zero_jet': dummy_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14]
    data = np.delete(data, dummy_list, 1)

    names = [
             "jet1_pt", "jet1_eta", "jet1_phi*",
             "jet2_pt", "jet2_eta", "jet2_phi*",
             "met",     "met_phi",
             "dijet_pt", "dijet_eta", "dijet_phi*", "dijet_m",
             "Z_pt", "Z_eta", "dummy", "costheta*"
            ]

    if region == 'two_jet': names += ["Higgs score", "VBF score", "m_mumu"]
    else: names += ["Higgs score", "m_mumu"]

    names = np.delete(names, dummy_list, 0)
    print names


    data = pd.DataFrame(data=data, columns=names)
    corr = data.corr()

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask, k=1)] = True

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    plt.figure(figsize=(18,18))
    sns.set(font_scale=2)

    ax = sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1.0, vmin=-1.0, xticklabels = 1, yticklabels = 1, square=True, linewidth=0.5)
    ax.set_title(title + '_' + region)

    #plt.show()
    plt.savefig('corr/' + title + '_' + region + ".png")


    return



def main():
    args=getArgs()
    

    sigs = ['ggF','VBF','VH','ttH']
    bkgs = ['data']

    region = args.region

    apply_weight(sigs, region, 'signals')
    apply_weight(bkgs, region, 'backgrounds')

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
