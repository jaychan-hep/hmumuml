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
from tqdm import tqdm
import logging
from pdb import set_trace
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

pd.options.mode.chained_assignment = None

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', nargs=2, default=['data/training_config.json', 'data/apply_config.json'], help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', default='skimmed_ntuples', help='directory of training inputs')
    parser.add_argument('-m', '--modelFolder', action='store', default='models', help='directory of BDT models')
    parser.add_argument('-o', '--outputFolder', action='store', default='outputs', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'all_jet'], default='zero_jet', help='Region to process')
    return parser.parse_args()

class ApplyXGBHandler(object):
    "Class for applying XGBoost"

    def __init__(self, configPath, region=''):

        print '==============================='
        print '  ApplyXGBHandler initialized'
        print '==============================='

        self._region = region
        self._inputFolder = ''
        self._inputTree = region if region else 'inclusive'
        self._modelFolder = ''
        self._outputFolder = ''
        self._chunksize = 500000
        self._category = []
        self._branches = []
        self._outbranches = []

        self.m_models = {}
        self.m_tsfs = {}

        self.train_variables = {}
        self.randomIndex = 'eventNumber'

        self.models = {}
        self.observables = []
        self.preselections = []

        self.readApplyConfig(configPath[1])
        self.readTrainConfig(configPath[0])
        self.arrangeBranches()
        self.arrangePreselections()

    def readApplyConfig(self, configPath):
        """Read configuration file formated in json to extract information to fill TemplateMaker variables."""
        try:
            member_variables = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("_") and not attr.startswith('m_')]

            stream = open(configPath, 'r')
            configs = json.loads(stream.read())

            # read from the common settings
            config = configs["common"]
            for member in config.keys():
                if member in member_variables:
                    setattr(self, member, config[member])

            # read from the region specific settings
            if self._region:
                config = configs[self._region]
                for member in config.keys():
                    if member in member_variables:
                        setattr(self, member, config[member])
                if '+preselections' in config.keys():
                    self.preselections += config['+preselections']
                if '+observables' in config.keys():
                    self.observables += config['+observables']

        except Exception as e:
            logging.error("Error reading apply configuration '{config}'".format(config=configPath))
            logging.error(e)

    def readTrainConfig(self, configPath):

        try:
            stream = open(configPath, 'r')
            configs = json.loads(stream.read())
   
            config = configs["common"]
            if 'randomIndex' in config.keys(): self.randomIndex = config['randomIndex']
 
            if self.models:
                for model in self.models:
    
                    # read from the common settings
                    config = configs["common"]
                    if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables'][:]
    
                    # read from the region specific settings
                    if model in configs.keys():
                        config = configs[model]
                        if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables'][:]
                        if '+train_variables' in config.keys(): self.train_variables[model] += config['+train_variables']

        except Exception as e:
            logging.error("Error reading training configuration '{config}'".format(config=configPath))
            logging.error(e)

    def arrangeBranches(self):

        self._branches = set()
        for model in self.models:
            self._branches = self._branches | set(self.train_variables[model])

        self._branches = self._branches | set([self.randomIndex]) | set([p.split()[0] for p in self.preselections]) | set(self.observables)
        self._branches = list(self._branches)

        for model in self.models:
            self.train_variables[model] = [x.replace('noexpand:', '') for x in self.train_variables[model]]
        self.preselections = [x.replace('noexpand:', '') for x in self.preselections]
        self.randomIndex = self.randomIndex.replace('noexpand:', '')

        self._outbranches = [branch for branch in self._branches if 'noexpand' not in branch]

    def arrangePreselections(self):

        if self.preselections:
            self.preselections = ['data.' + p for p in self.preselections]

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setModelFolder(self, modelFolder):
        self._modelFolder = modelFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def preselect(self, data):

        for p in self.preselections:
            data = data[eval(p)]

        return data

    def loadModels(self):

        if self.models:
            for model in self.models:
                print 'XGB INFO: Loading BDT model: ', model
                self.m_models[model] = []
                for i in range(4):
                    bst = xgb.Booster()
                    bst.load_model('%s/%s_%d.h5'%(self._modelFolder, model, i))
                    self.m_models[model].append(bst)
                    del bst

    def loadTransformer(self):
        
        if self.models:
            for model in self.models:
                print 'XGB INFO: Loading score transformer for model: ', model
                self.m_tsfs[model] = []
                for i in range(4):
                    tsf = pickle.load(open('%s/tsf_%s_%d.pkl'%(self._modelFolder, model, i), "rb" ) )
                    self.m_tsfs[model].append(tsf)

    def applyBDT(self, category, scale=1):

        outputContainer = self._outputFolder + '/' + self._region
        output_path = outputContainer + '/%s.root' % category
        if not os.path.isdir(outputContainer): os.makedirs(outputContainer)
        if os.path.isfile(output_path): os.remove(output_path)

        f_list = []
        cat_folder = self._inputFolder + '/' + category
        for f in os.listdir(cat_folder):
            if f.endswith('.root'): f_list.append(cat_folder + '/' + f)

        print '-------------------------------------------------'
        for f in f_list: print 'XGB INFO: Including sample: ', f

        #TODO put this to the config
        for data in tqdm(read_root(sorted(f_list), key=self._inputTree, columns=self._branches, chunksize=self._chunksize), ncols=100, desc='XGB INFO: Applying BDTs to %s samples' % category):
            data = self.preselect(data)

            out_data = pd.DataFrame()

            for i in range(4):
                data_s = data[data[self.randomIndex]%4 == i]
                data_o = data_s[self._outbranches]

                for model in self.train_variables.keys():
                    x_Events = data_s[self.train_variables[model]]
                    dEvents = xgb.DMatrix(x_Events)
                    scores = self.m_models[model][i].predict(dEvents)
                    scores_t = self.m_tsfs[model][i].transform(scores.reshape(-1,1)).reshape(-1)
                
                    xgb_basename = self.models[model]
                    data_o[xgb_basename] = scores
                    data_o[xgb_basename+'_t'] = scores_t

                out_data = pd.concat([out_data, data_o], ignore_index=True, sort=False)

            out_data.to_root(output_path, key='test', mode='a', index=False)

            del out_data

def main():

    args=getArgs()
    
    configPath = args.config
    xgb = ApplyXGBHandler(configPath, args.region)

    xgb.setInputFolder(args.inputFolder)
    xgb.setModelFolder(args.modelFolder)
    xgb.setOutputFolder(args.outputFolder)

    xgb.loadModels()
    xgb.loadTransformer()

    with open('data/inputs_config.json') as f:
        config = json.load(f)
    sample_list = config['sample_list']

    for category in sample_list:
        xgb.applyBDT(category)

    return

if __name__ == '__main__':
    main()
