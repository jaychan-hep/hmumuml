#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.21.2019
#
#
#
#
import os
import math
from argparse import ArgumentParser
from ROOT import *
#import numpy as np
#import time
import pandas as pd
from root_pandas import *
from progressbar import ProgressBar

def getArgs():
    parser = ArgumentParser(description="Skim the input ntuples for Hmumu XGBoost analysis.")
    parser.add_argument('-i', '--input', action='store', default='inputs', help='Path to the input ntuple')
    parser.add_argument('-o', '--output', action='store', default='outputs', help='Path to the output ntuple')
    parser.add_argument('--chunksize', type=int, default=500000, help='size to process at a time') 
    return  parser.parse_args()

def compute_Y_jj(x):

    if x.Jets_jetMultip < 2:
        return -9999
    else:
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiE4D<float>")(x.Jets_PT_Lead, x.Jets_Eta_Lead, x.Jets_Phi_Lead, x.Jets_E_Lead)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiE4D<float>")(x.Jets_PT_Sub, x.Jets_Eta_Sub, x.Jets_Phi_Sub, x.Jets_E_Sub)
        j1j2 = Jets_1+Jets_2
        return j1j2.Rapidity()

def compute_Delta_Phi(x, var, min_jet=0):

    if x.Jets_jetMultip < min_jet: return -9999
    return TVector2.Phi_mpi_pi(x[var] - x.Z_Phi_FSR)

def preselect(data):

    data = data[data.Muons_Minv_MuMu_Fsr >= 110]
    #data = data[data.Jets_jetMultip >= 2]
    data = data[data.FinalSelection == True]
    data = data[data.SampleOverlapWeight == True]
    data = data[data.EventWeight_MCCleaning_5 != 0]

    return data

def decorate(data):

    if data.shape[0] == 0: return data

    data['weight'] = data.GlobalWeight * data.SampleOverlapWeight
    data['Jets_Y_jj'] = data.apply(lambda x: compute_Y_jj(x), axis=1)
    data['DeltaPhi_mumumu1'] = data.apply(lambda x: compute_Delta_Phi(x, 'Muons_Phi_Lead'), axis=1)
    data['DeltaPhi_mumumu2'] = data.apply(lambda x: compute_Delta_Phi(x, 'Muons_Phi_Sub'), axis=1)
    data['DeltaPhi_mumumj1'] = data.apply(lambda x: compute_Delta_Phi(x, 'Jets_Phi_Lead', min_jet=1), axis=1)
    data['DeltaPhi_mumumj2'] = data.apply(lambda x: compute_Delta_Phi(x, 'Jets_Phi_Sub', min_jet=2), axis=1)
    data['DeltaPhi_mumujj'] = data.apply(lambda x: compute_Delta_Phi(x, 'Jets_Phi_jj', min_jet=2), axis=1)
    data['DeltaPhi_mumuMET'] = data.apply(lambda x: compute_Delta_Phi(x, 'Event_MET_Phi'), axis=1)

    return data
    

def main():

    args = getArgs()

    variables = ['EventInfo_EventNumber', 'FinalSelection', 
                 'Muons_*', 'Z_*', 'Jets_*', 'Event_MET', 'Event_MET_Phi',
                 'GlobalWeight', 'SampleOverlapWeight', 'EventWeight_MCCleaning_5']

    if os.path.isfile(args.output): os.remove(args.output)

    initial_events = 0
    final_events = 0

    pbar = ProgressBar()
    for data in pbar(read_root(args.input, key='DiMuonNtuple', columns=variables, chunksize=args.chunksize)):

        initial_events += data.shape[0]
        #data = preprocess(data)
        data = preselect(data) #TODO add cutflow
        data = decorate(data)
        final_events += data.shape[0]
        data.rename(columns={'Muons_Minv_MuMu_Fsr': 'm_mumu', 'EventInfo_EventNumber': 'eventNumber', 'Jets_jetMultip': 'n_j'}, inplace=True)
        data.to_root(args.output, key='Skimmed_Hmumu', mode='a', index=False)

    meta_data = pd.DataFrame({'initial_events': [initial_events], 'final_events': [final_events]})
    meta_data.to_root(args.output, key='MetaData', mode='a', index=False)

if __name__ == '__main__':
    main()
