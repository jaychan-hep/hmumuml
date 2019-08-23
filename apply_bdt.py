#!/usr/bin/env python
import os
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd
from root_pandas import *
import pickle
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import xgboost as xgb
from progressbar import ProgressBar
import logging
from pdb import set_trace
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', default='data/training_config.json', help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', help='directory of training inputs')
    parser.add_argument('-m', '--modelFolder', action='store', help='directory of BDT models')
    parser.add_argument('-o', '--outputFolder', action='store', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'all_jet'], default='zero_jet', help='Region to process')
    return parser.parse_args()

class ApplyXGBHandler(object):
    "Class for applying XGBoost"

    def __init__(self, configPath, region=''):

        self._region = region
        self._inputFolder = 'skimmed_ntuples'
        self._inputTree = 'Skimmed_Hmumu'
        self._modelFolder = 'models'
        self._outputFolder = 'outputs'
        self._outputContainer = self._outputFolder + '/' + self._region
        self._chunksize = 500000
        self._randomIndex = 'eventNumber'
        self._weight = 'weight'
        self._category = []
        self.configPath = configPath

        self.m_models = {}
        self.m_tsfs = {}

        self.train_variables = {}

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def preselect(self, data):

        #TODO put this to the config
        if self._region == 'zero_jet': data = data[data.n_j == 0]
        elif self._region == 'one_jet': data = data[data.n_j == 1]
        elif self._region in ['two_jet', 'VBF']: data = data[data.n_j >= 2]

        return data

    def loadModels(self, models=[]):

        stream = open(self.configPath, 'r')
        configs = json.loads(stream.read())

        if models:
            for model in models:
                self.m_models[model] = []
                for i in range(4):
                    bst = xgb.Booster()
                    bst.load_model('%s/%s_%d.h5'%(self._modelFolder, model, i))
                    self.m_models[model].append(bst)
                    del bst

                # read from the common settings
                config = configs["common"]
                if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables']
                    
                # read from the region specific settings
                if model in configs.keys():
                    config = configs[model]
                    if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables']
                    if '+train_variables' in config.keys(): self.train_variables[model] += config['+train_variables']


    def loadTransformer(self, models=[]):
        
        if models:
            for model in models:
                self.m_tsfs[model] = []
                for i in range(4):
                    tsf = pickle.load(open('%s/tsf_%s_%d.pkl'%(self._modelFolder, model, i), "rb" ) )
                    self.m_tsfs[model].append(tsf)

    def applyBDT(self, category):

        output_path = self._outputContainer + '/%s.root' % category
        if not os.path.isdir(self._outputContainer): os.makedirs(self._outputContainer)
        if os.path.isfile(output_path): os.remove(output_path)

        branches = set()
        for model in self.train_variables.keys():
            branches = branches | set(self.train_variables[model])

        branches = branches | set(['m_mumu', 'n_j']) | set([self._randomIndex]) | set([self._weight])
        branches = list(branches)

        f_list = []
        cat_folder = self._inputFolder + '/' + category
        for f in os.listdir(cat_folder):
            if f.endswith('.root'): f_list.append(cat_folder + '/' + f)

        print '-------------------------------------------------'
        print 'XGB INFO: Applying BDTs to %s samples...' % category
        for f in f_list: print 'XGB INFO: Including sample: ', f
        pbar = ProgressBar()
        #TODO put this to the config
        for data in pbar(read_root(sorted(f_list), key=self._inputTree, columns=branches, chunksize=self._chunksize)):
            data = self.preselect(data)
            for model in self.train_variables.keys():
                index_Events = data[self._randomIndex]
                x_Events = data[self.train_variables[model]]
                dEvents = xgb.DMatrix(x_Events.values)
                for i in range(4):
                    scores = self.m_models[model][i].predict(dEvents)
                    scores_t = self.m_tsfs[model][i].transform(scores.reshape(-1,1)).reshape(-1)
                    
                    index_Events = pd.concat([index_Events, pd.DataFrame(index=index_Events.index,data = {'xgb_%d' % i: scores,'tsf_%d' % i: scores_t})], axis=1)

                xgb_basename = 'bdt_score' if model != 'VBF' else 'bdt_score_VBF'
                score_Events = index_Events.apply(lambda x: self.getScore(x) , axis=1, result_type='expand')
                score_Events.columns = [xgb_basename, xgb_basename+'_t', xgb_basename+'_val', xgb_basename+'_val_t', xgb_basename+'_trainA', xgb_basename+'_trainA_t', xgb_basename+'_trainB', xgb_basename+'_trainB_t']
                data = pd.concat([data, score_Events], axis=1)
                data.to_root(output_path, key='test', mode='a', index=False)


    def getScore(self, x):

        tag = x[self._randomIndex]%4
        tag_val = (tag - 1)%4
        tag_trainA = (tag - 2)%4
        tag_trainB = (tag - 3)%4

        bdt_score = x['xgb_%d' % tag]
        bdt_score_t = x['tsf_%d' % tag]
        bdt_score_val = x['xgb_%d' % tag_val]
        bdt_score_val_t = x['tsf_%d' % tag_val]
        bdt_score_trainA = x['xgb_%d' % tag_trainA]
        bdt_score_trainA_t = x['tsf_%d' % tag_trainA]
        bdt_score_trainB = x['xgb_%d' % tag_trainB]
        bdt_score_trainB_t = x['tsf_%d' % tag_trainB]

        return bdt_score, bdt_score_t, bdt_score_val, bdt_score_val_t, bdt_score_trainA, bdt_score_trainA_t, bdt_score_trainB, bdt_score_trainB_t 

def main():

    args=getArgs()
    
    configPath = args.config
    xgb = ApplyXGBHandler(configPath, args.region)

    models = [args.region]
    if args.region == 'two_jet': models.append('VBF')
    xgb.loadModels(models)
    xgb.loadTransformer(models)

#    with open('data/inputs_config.json') as f:
#        config = json.load(f)
#    sample_list = config['sample_list']

    categories = []
    categories += ['ggF','VBF','VH','ttH']
    categories += ['data_sid']
    categories += ['Z', 'ttbar', 'diboson', 'stop']

    for category in categories:
        xgb.applyBDT(category)

    return

if __name__ == '__main__':
    main()
