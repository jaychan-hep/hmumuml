#!/usr/bin/env python
import os
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd
from root_pandas import *
import pickle
import random
from sklearn.metrics import roc_curve, auc, confusion_matrix, roc_auc_score
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import xgboost as xgb
from tabulate import tabulate
#from bayes_opt import BayesianOptimization
import matplotlib.pyplot as plt
from tqdm import tqdm
#import logging
from pdb import set_trace
#logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
#logging.getLogger('matplotlib.font_manager').disabled = True
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError + 1

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import AveragePooling2D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import concatenate
from tensorflow.keras.layers.experimental.preprocessing import Normalization
import tensorflow.keras as keras

from sklearn.preprocessing import StandardScaler
from pickle import dump, load

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', default='data/training_config_NN.json', help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', help='directory of training inputs')
    parser.add_argument('-o', '--outputFolder', action='store', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'VBF', 'all_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-f', '--fold', action='store', type=int, nargs='+', choices=[0, 1, 2, 3], default=[0, 1, 2, 3], help='specify the fold for training')
    parser.add_argument('-p', '--params', action='store', type=json.loads, default=None, help='json string.') #type=json.loads
    parser.add_argument('--save', action='store_true', help='Save model weights to HDF5 file')
    parser.add_argument('--roc', action='store_true', help='Plot ROC')
    parser.add_argument('--skopt', action='store_true', default=False, help='Run hyperparameter tuning using skopt')
    parser.add_argument('--skopt-plot', action='store_true', default=False, help='Plot skopt results')
    return parser.parse_args()

def fixRandomSeed(seed_value):

    os.environ['PYTHONHASHSEED']=str(seed_value)
    random.seed(seed_value)
    np.random.seed(seed_value)
    tf.compat.v1.set_random_seed(seed_value)
    tf.random.set_seed(seed_value)

class XGBoostHandler(object):
    "Class for running XGBoost"

    def __init__(self, configPath, region=''):

        print("""
========================================================================
|| **     *     ** ###### $$        ######    ****    #       # %%%%% ||
||  **   ***   **  ##     $$       ##       ***  ***  ##     ## %     ||
||   ** ** ** **   ###### $$      ##       **      ** # #   # # %%%%% ||
||    ***   ***    ##     $$       ##       ***  ***  #  # #  # %     || 
||     *     *     ###### $$$$$$$   ######    ****    #   #   # %%%%% || 
||                                                                    ||
||    $$$$$$$$$$    ####         XX      XX     GGGGGG     BBBBB      ||
||        $$      ###  ###         XX  XX     GGG          B    BB    ||
||        $$     ##      ##          XX      GG     GGGG   BBBBB      ||
||        $$      ###  ###         XX  XX     GGG     GG   B    BB    ||
||        $$        ####         XX      XX     GGGGGG     BBBBB      ||
========================================================================
              """)


        self._region = region
        self._inputFolder = 'skimmed_ntuples'
        self._outputFolder = 'models'
        self._chunksize = 500000
        self._branches = []

        self.m_data_sig = pd.DataFrame()
        self.m_data_bkg = pd.DataFrame()
        self.m_train_wt = {}
        self.m_val_wt = {}
        self.m_test_wt = {}
        self.m_x_train = {}
        self.m_x_val = {}
        self.m_x_test = {}
        self.m_x_test_sig = {}
        self.m_x_test_bkg = {}
        self.m_y_train = {}
        self.m_y_val = {}
        self.m_y_test = {}
        self.m_model = {}
        self.m_metrics = {}
        self.m_score_train = {}
        self.m_score_val = {}
        self.m_score_test = {}
        self.m_score_test_sig = {}
        self.m_score_test_bkg = {}
        self.m_tsf = {}

        self.inputTree = 'inclusive'
        self.train_signal = []
        self.train_background = []
        self.train_variables = []
        self.preselections = []
        self.signal_preselections = []
        self.background_preselections = []
        self.randomIndex = 'eventNumber'
        self.weight = 'weight'
        self.params = [{'eval_metric': ['auc', 'logloss']}]
        self.early_stopping_rounds = 10
        self.numRound = 10000
        self.SF = 1.

        self.readConfig(configPath)
        self.checkConfig()

    def readConfig(self, configPath):
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
                if '+train_variables' in config.keys():
                    self.train_variables += config['+train_variables']
                if '+preselections' in config.keys():
                    self.preselections += config['+preselections']
                if '+signal_preselections' in config.keys():
                    self.signal_preselections += config['+signal_preselections']
                if '+background_preselections' in config.keys():
                    self.background_preselections += config['+background_preselections']

            self._branches = list( set(self.train_variables) | set([p.split()[0] for p in self.preselections]) | set([p.split()[0] for p in self.signal_preselections]) | set([p.split()[0] for p in self.background_preselections]) | set([self.randomIndex, self.weight]))

            self.train_variables = [x.replace('noexpand:', '') for x in self.train_variables]
            self.preselections = [x.replace('noexpand:', '') for x in self.preselections]
            self.signal_preselections = [x.replace('noexpand:', '') for x in self.signal_preselections]
            self.background_preselections = [x.replace('noexpand:', '') for x in self.background_preselections]
            self.randomIndex = self.randomIndex.replace('noexpand:', '')
            self.weight = self.weight.replace('noexpand:', '')

            if self.preselections:
                self.preselections = ['data.' + p for p in self.preselections]
            if self.signal_preselections:
                self.signal_preselections = ['data.' + p for p in self.signal_preselections]
            if self.background_preselections:
                self.background_preselections = ['data.' + p for p in self.background_preselections]

        except Exception as e:
            logging.error("Error reading configuration '{config}'".format(config=configPath))
            logging.error(e)

    def checkConfig(self):
        if not self.train_signal: print('ERROR: no training signal!!')
        if not self.train_background: print('ERROR: no training background!!')
        if not self.train_variables: print('ERROR: no training variables!!')

    def setParams(self, params, fold=-1):

        print('XGB INFO: setting hyperparameters...')
        if fold == -1:
            self.params = [{'eval_metric': ['auc', 'logloss']}]
            fold = 0
        for key in params:
            self.params[fold][key] = params[key]

    def set_early_stopping_rounds(self, rounds):
        self.early_stopping_rounds = rounds

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def preselect(self, data, sample=''):

        if sample == 'signal':
            for p in self.signal_preselections:
                data = data[eval(p)]
        elif sample == 'background':
            for p in self.background_preselections:
                data = data[eval(p)]
        for p in self.preselections:
            data = data[eval(p)]

        return data

    def readData(self):

        sig_list, bkg_list = [], []
        for sig_cat in self.train_signal:
            sig_cat_folder = self._inputFolder + '/' + sig_cat
            for sig in os.listdir(sig_cat_folder):
                if sig.endswith('.root'): sig_list.append(sig_cat_folder + '/' + sig)
        for bkg_cat in self.train_background:
            bkg_cat_folder = self._inputFolder + '/' + bkg_cat
            for bkg in os.listdir(bkg_cat_folder):
                if bkg.endswith('.root'): bkg_list.append(bkg_cat_folder + '/' + bkg)

        print('-------------------------------------------------')
        for sig in sig_list: print('XGB INFO: Adding signal sample: ', sig)
        #TODO put this to the config
        for data in tqdm(read_root(sorted(sig_list), key=self.inputTree, columns=self._branches, chunksize=self._chunksize), desc='XGB INFO: Loading training signals', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
            data = self.preselect(data, 'signal')
            self.m_data_sig = self.m_data_sig.append(data, ignore_index=True)

        print('----------------------------------------------------------')
        for bkg in bkg_list: print('XGB INFO: Adding background sample: ', bkg)
        #TODO put this to the config
        for data in tqdm(read_root(sorted(bkg_list), key=self.inputTree, columns=self._branches, chunksize=self._chunksize), desc='XGB INFO: Loading training backgrounds', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):
            data = self.preselect(data, 'background')
            self.m_data_bkg = self.m_data_bkg.append(data, ignore_index=True)

    def prepareData(self, fold=0):

        # training, validation, test split
        print('----------------------------------------------------------')
        print('XGB INFO: Splitting samples to training, validation, and test...')
        test_sig = self.m_data_sig[self.m_data_sig[self.randomIndex]%4 == fold]
        test_bkg = self.m_data_bkg[self.m_data_bkg[self.randomIndex]%4 == fold]
        val_sig = self.m_data_sig[(self.m_data_sig[self.randomIndex]-1)%4 == fold]
        val_bkg = self.m_data_bkg[(self.m_data_bkg[self.randomIndex]-1)%4 == fold]
        train_sig = self.m_data_sig[((self.m_data_sig[self.randomIndex]-2)%4 == fold) | ((self.m_data_sig[self.randomIndex]-3)%4 == fold)]
        train_bkg = self.m_data_bkg[((self.m_data_bkg[self.randomIndex]-2)%4 == fold) | ((self.m_data_bkg[self.randomIndex]-3)%4 == fold)]

        headers = ['Sample', 'Total', 'Training', 'Validation']
        sample_size_table = [
            ['Signal'    , len(train_sig)+len(val_sig), len(train_sig), len(val_sig)],
            ['Background', len(train_bkg)+len(val_bkg), len(train_bkg), len(val_bkg)],
        ]

        print(tabulate(sample_size_table, headers=headers, tablefmt='simple'))

        # setup the training data set
        print('XGB INFO: Setting the training arrays...')
        x_test_sig = test_sig[self.train_variables]
        x_test_bkg = test_bkg[self.train_variables]
        x_val_sig = val_sig[self.train_variables]
        x_val_bkg = val_bkg[self.train_variables]
        x_train_sig = train_sig[self.train_variables]
        x_train_bkg = train_bkg[self.train_variables]

        x_train = pd.concat([x_train_sig, x_train_bkg])
        x_val = pd.concat([x_val_sig, x_val_bkg])
        x_test = pd.concat([x_test_sig, x_test_bkg])
        
        self.m_x_train[fold] = x_train.to_numpy()
        self.m_x_val[fold] = x_val.to_numpy()
        self.m_x_test[fold] = x_test.to_numpy()
   
        self.m_x_test_sig[fold] = x_test_sig.to_numpy()
        self.m_x_test_bkg[fold] = x_test_bkg.to_numpy()

        # setup the weights
        print('XGB INFO: Setting the event weights...')

        if self.SF == -1: self.SF = 1.*len(train_sig_wt)/len(train_bkg_wt)

        train_sig_wt = train_sig[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_sig[self.weight].mean() * (1.+self.SF) * train_sig.shape[0] )
        train_bkg_wt = train_bkg[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_bkg[self.weight].mean() * (1.+self.SF) * train_bkg.shape[0] )
        val_sig_wt = val_sig[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_sig[self.weight].mean() * (1.+self.SF) * train_sig.shape[0] )
        val_bkg_wt = val_bkg[[self.weight]] * ( train_sig.shape[0] + train_bkg.shape[0] ) * self.SF / ( train_bkg[self.weight].mean() * (1.+self.SF) * train_bkg.shape[0] )
        test_sig_wt = test_sig[[self.weight]]
        test_bkg_wt = test_bkg[[self.weight]]

        self.m_train_wt[fold] = pd.concat([train_sig_wt, train_bkg_wt]).to_numpy().reshape((-1))
        self.m_val_wt[fold] = pd.concat([val_sig_wt, val_bkg_wt]).to_numpy().reshape((-1))
        self.m_test_wt[fold] = pd.concat([test_sig_wt, test_bkg_wt]).to_numpy().reshape((-1))

        # setup the truth labels
        print('XGB INFO: Signal labeled as one; background labeled as zero.')
        self.m_y_train[fold] = np.concatenate((np.ones(len(train_sig), dtype=np.uint8), np.zeros(len(train_bkg), dtype=np.uint8)))
        self.m_y_val[fold]   = np.concatenate((np.ones(len(val_sig)  , dtype=np.uint8), np.zeros(len(val_bkg)  , dtype=np.uint8)))
        self.m_y_test[fold]  = np.concatenate((np.ones(len(test_sig)  , dtype=np.uint8), np.zeros(len(test_bkg)  , dtype=np.uint8)))

    def AddMetric(self, metric):

        if metric == "auc":
            self.m_metrics["auc"] = keras.metrics.AUC(name="auc")

    def buildModel(self, fold, param=None):

        if not param:
            print('XGB INFO: Setting the hyperparameters...')
            param = self.params[0 if len(self.params) == 1 else fold]
        print('param: ', param)

        print('XGB INFO: Building neural network architecture...')
        input_shape = self.m_x_train[fold].shape[1]
        print(f"Input shape: {input_shape}")

        if not "auc" in self.m_metrics:
            self.AddMetric("auc")

        normalizer = Normalization()
        normalizer.adapt(self.m_x_train[fold])

        self.m_model[fold] = Sequential()
        self.m_model[fold].add(Input(shape=(input_shape,)))
        self.m_model[fold].add(normalizer)
        self.m_model[fold].add(Dense(units=64, activation='relu'))
        self.m_model[fold].add(Dense(units=32, activation='relu'))
        self.m_model[fold].add(Dropout(0.1))
        self.m_model[fold].add(Dense(units=16, activation='relu'))
        self.m_model[fold].add(Dense(units=8, activation='relu'))
        self.m_model[fold].add(Dense(units=1, activation='sigmoid'))

        self.m_model[fold].compile(loss=keras.losses.BinaryCrossentropy(),
              optimizer=tf.keras.optimizers.Adam(param.get("lr", 0.001)),
              # optimizer=keras.optimizers.SGD(learning_rate=param.get("lr", 0.001), momentum=param.get("momentum", 0.9), nesterov=True),
              weighted_metrics=[self.m_metrics["auc"]])
      
        self.m_model[fold].summary()

    def trainModel(self, fold, param=None):

        if not param:
            print('XGB INFO: Setting the hyperparameters...')
            param = self.params[0 if len(self.params) == 1 else fold]
        print('param: ', param)

        # finally start training!!!
        print('XGB INFO: Start training!!!')
        callback = keras.callbacks.EarlyStopping(monitor='val_auc', patience=self.early_stopping_rounds, restore_best_weights=True, mode='max')
        try:
            self.m_model[fold].fit(self.m_x_train[fold],
                self.m_y_train[fold],
                epochs=self.numRound,
                batch_size=param.get("Bs", 64),
                sample_weight=self.m_train_wt[fold],
                validation_data=(self.m_x_val[fold], self.m_y_val[fold], self.m_val_wt[fold]),
                callbacks=[callback]
            )
        except KeyboardInterrupt:
            print('Finishing on SIGINT.')

    def testModel(self, fold):
        # get scores
        print('XGB INFO: Computing scores for different sample sets...')
        self.m_score_val[fold]   = self.m_model[fold].predict(self.m_x_val[fold])
        self.m_score_test[fold]  = self.m_model[fold].predict(self.m_x_test[fold])
        self.m_score_train[fold] = self.m_model[fold].predict(self.m_x_train[fold])
        self.m_score_test_sig[fold]  = self.m_model[fold].predict(self.m_x_test_sig[fold])
        self.m_score_test_bkg[fold]  = self.m_model[fold].predict(self.m_x_test_bkg[fold])

    def plotScore(self, fold=0, sample_set='val'):

        if sample_set == 'train':
            plt.hist(self.m_score_train[fold], bins='auto')
        elif sample_set == 'val':
            plt.hist(self.m_score_val[fold], bins='auto')
        elif sample_set == 'test':
            plt.hist(self.m_score_test[fold], bins='auto')
        plt.title('Score distribution %s' % sample_set)
        plt.show()
        if sample_set == 'test':
            plt.hist(self.m_score_test_sig[fold], bins='auto')
            plt.title('Score distribution test signal')
            plt.show()
            plt.hist(self.m_score_test_bkg[fold], bins='auto')
            plt.title('Score distribution test background')
            plt.show()

    def plotFeaturesImportance(self, fold=0, save=True, show=False, type='weight'):
        """Plot feature importance. Type can be 'weight', 'gain' or 'cover'"""

        xgb.plot_importance(booster=self.m_bst[fold], importance_type='weight')
        plt.tight_layout()

        if save:
            # create output directory
            os.makedirs('plots/feature_importance', exist_ok=True)
            # save figure
            plt.savefig('plots/feature_importance/%s_%d.pdf' % (self._region, fold))

        if show: plt.show()

        plt.clf()

        return

    def plotROC(self, fold=0, save=True, show=False):
    
        fpr_train, tpr_train, _ = roc_curve(self.m_y_train[fold], self.m_score_train[fold], sample_weight=self.m_train_wt[fold])
        tpr_train, fpr_train = np.array(list(zip(*sorted(zip(tpr_train, fpr_train)))))
        roc_auc_train = 1 - auc(tpr_train, fpr_train)
        fpr_val, tpr_val, _ = roc_curve(self.m_y_val[fold], self.m_score_val[fold], sample_weight=self.m_val_wt[fold])
        tpr_val, fpr_val = np.array(list(zip(*sorted(zip(tpr_val, fpr_val)))))
        roc_auc_val = 1 - auc(tpr_val, fpr_val)
        fpr_test, tpr_test, _ = roc_curve(self.m_y_test[fold], self.m_score_test[fold], sample_weight=self.m_test_wt[fold])
        tpr_test, fpr_test = np.array(list(zip(*sorted(zip(tpr_test, fpr_test)))))
        roc_auc_test = 1 - auc(tpr_test, fpr_test)

        fnr_train = 1.0 - fpr_train
        fnr_val = 1.0 - fpr_val
        fnr_test = 1.0 - fpr_test

        plt.grid(color='gray', linestyle='--', linewidth=1)
        plt.plot(tpr_train, fnr_train, label='Train set, area = %0.6f' % roc_auc_train, color='black', linestyle='dotted')
        plt.plot(tpr_val, fnr_val, label='Val set, area = %0.6f' % roc_auc_val, color='blue', linestyle='dashdot')
        plt.plot(tpr_test, fnr_test, label='Test set, area = %0.6f' % roc_auc_test, color='red', linestyle='dashed')
        plt.plot([0, 1], [1, 0], linestyle='--', color='black', label='Luck')
        plt.xlabel('Signal acceptance')
        plt.ylabel('Background rejection')
        plt.title('Receiver operating characteristic')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xticks(np.arange(0, 1, 0.1))
        plt.yticks(np.arange(0, 1, 0.1))
        plt.legend(loc='lower left', framealpha=1.0)
        plt.tight_layout()

        if save:
            # create output directory
            os.makedirs('plots/roc_curve', exist_ok=True)
            # save figure
            plt.savefig('plots/roc_curve/NN_%s_%d.pdf' % (self._region, fold))

        if show: plt.show()

        plt.clf()

        return

    def getAUC(self, fold=0, sample_set='val'):

        if sample_set == 'train':
            fpr, tpr, _ = roc_curve(self.m_y_train[fold], self.m_score_train[fold], sample_weight=self.m_train_wt[fold])
        elif sample_set == 'val':
            fpr, tpr, _ = roc_curve(self.m_y_val[fold], self.m_score_val[fold], sample_weight=self.m_val_wt[fold])
        elif sample_set == 'test':
            fpr, tpr, _ = roc_curve(self.m_y_test[fold], self.m_score_test[fold], sample_weight=self.m_test_wt[fold])

        tpr, fpr = np.array(list(zip(*sorted(zip(tpr, fpr)))))
        roc_auc = 1 - auc(tpr, fpr)

        return self.params[0 if len(self.params) == 1 else fold], roc_auc

    def transformScore(self, fold=0, sample='sig'):

        print(f'XGB INFO: transforming scores based on {sample}')
        # transform the scores
        self.m_tsf[fold] = QuantileTransformer(n_quantiles=1000, output_distribution='uniform', subsample=1000000000, random_state=0)
        #plt.hist(score_test_sig, bins='auto')
        #plt.show()
        if sample == 'sig': self.m_tsf[fold].fit(self.m_score_test_sig[fold].reshape(-1, 1))
        elif sample == 'bkg': self.m_tsf[fold].fit(self.m_score_test_bkg[fold].reshape(-1, 1))
        #score_test_sig_t=tsf.transform(self.m_score_test_sig[fold].reshape(-1, 1)).reshape(-1)
        #plt.hist(score_test_sig_t, bins='auto')
        #plt.show()

    def save(self, fold=0):

        # create output directory
        os.makedirs(self._outputFolder, exist_ok=True)

        # save BDT model
        if fold in self.m_model.keys(): self.m_model[fold].save('%s/NN_%s_%d.tf' % (self._outputFolder, self._region, fold))

        # save score transformer
        if fold in self.m_tsf.keys(): 
            with open('%s/tsf_%s_%d.pkl' %(self._outputFolder, self._region, fold), 'wb') as f:
                pickle.dump(self.m_tsf[fold], f, -1)

    def skoptHP(self, fold):

        import skopt
        from skopt import dump
        from skopt import Optimizer
        from skopt import gbrt_minimize, gp_minimize
        from skopt.utils import use_named_args
        from skopt.space import Real, Categorical, Integer

        print ('default param: ', self._region, fold, self.params[fold])

        if self._region == "zero_jet":
            search_space = [
                Real(0.0001, 0.001, name='lr'),#, prior='log-uniform'),
                # Real(0., 1., name='momentum')
            ]
        elif self._region == "one_jet":
            search_space = [
                Real(0.0001, 0.001, name='lr'),#, prior='log-uniform'),
                # Real(0., 1., name='momentum')
            ]
        elif self._region == "two_jet":
            search_space = [
                Real(0.0001, 0.001, name='lr'),#, prior='log-uniform'),
                # Real(0., 1., name='momentum')
            ]
        elif self._region == "VBF":
            search_space = [
                Real(0.0001, 0.001, name='lr'),#, prior='log-uniform'),
                # Real(0., 1., name='momentum')
            ]

        def objective(param):

            self.setParams(param, fold)
            self.setParams({'eval_metric': ['auc']}, fold)
            self.buildModel(fold, self.params[fold])
            self.trainModel(fold, self.params[fold])
            self.testModel(fold)
            auc = self.getAUC(fold)[-1]
            return auc

        search_names = [var.name for var in search_space]

        opt = Optimizer(search_space, # TODO: Add noise
                    n_initial_points=10,
                    acq_optimizer_kwargs={'n_jobs':4})

        n_calls = 100
        exp_dir = 'models/skopt/'
        os.makedirs(exp_dir, exist_ok=True)

        best_AUC = 0
        for i in range(n_calls):

            sampled_point = opt.ask()
            param = dict(zip(search_names, sampled_point))
            # param = {key.decode(): val for key, val in param.items()}
            print('Point:', i+1, param)
            f_val = objective(param)
            print ('AUC:', f_val)
            if f_val > best_AUC:
                best_AUC = f_val
                self.plotROC(fold)
                self.transformScore(fold)
                self.save(fold)
            opt_result = opt.tell(sampled_point, -f_val)

            with open(exp_dir + 'optimizer_region_'+self._region+'_fold'+str(fold) + '.pkl', 'wb') as fp:
                
                dump(opt, fp)

            with open(exp_dir + 'results_region_'+self._region+'_fold'+str(fold) + '.pkl', 'wb') as fp:
                # Delete objective function before saving. To be used when results have been loaded,
                # the objective function must then be imported from this script.

                if opt_result.specs is not None:

                    if 'func' in opt_result.specs['args']:
                        res_without_func = copy.deepcopy(opt_result)
                        del res_without_func.specs['args']['func']
                        dump(res_without_func, fp)
                    else:
                        dump(opt_result, fp)

                else:

                    dump(opt_result, fp)

        with open(exp_dir + 'results_region_'+self._region+'_fold'+str(fold) + '.json', 'w') as fp:

            json.dump(dict(zip(search_names, opt_result.x)), fp)
            print('Save to ', exp_dir + 'results_region_'+self._region+'_fold'+str(fold) + '.json')


    def skoptPlot(self, fold):

        from skopt import load
        from skopt.plots import plot_objective, plot_evaluations, plot_convergence

        plt.switch_backend('agg') # For outputting plots when running on a server
        plt.ioff() # Turn off interactive mode, so that figures are only saved, not displayed
        plt.style.use('seaborn-paper')
        font_size = 9
        label_font_size = 8

        plt.rcParams.update({'font.size'        : font_size,
                 'axes.titlesize'   : font_size,
                 'axes.labelsize'   : font_size,
                 'xtick.labelsize'  : label_font_size,
                 'ytick.labelsize'  : label_font_size,
                 'figure.figsize'   : (5,4)}) # w,h

        exp_dir = 'models/skopt/'
        fig_dir = 'plots/skopt/figures/region_'+self._region+'_fold'+str(fold)+'/'

        os.makedirs(fig_dir, exist_ok=True)

        # Plots expect default rcParams
        plt.rcdefaults()

        # Load results
        res_loaded = load(exp_dir + 'results_region_'+self._region+'_fold'+str(fold) + '.pkl')
        print(res_loaded)

        # Plot evaluations
        fig,ax = plt.subplots()
        ax = plot_evaluations(res_loaded)
        plt.tight_layout()
        plt.savefig(fig_dir + 'evaluations.pdf')

        # Plot objective
        fig,ax = plt.subplots()
        ax = plot_objective(res_loaded)
        plt.tight_layout()
        plt.savefig(fig_dir + 'objective.pdf')

        # Plot convergence
        fig,ax = plt.subplots()
        ax = plot_convergence(res_loaded)
        plt.tight_layout()
        plt.savefig(fig_dir + 'convergence.pdf')

        ## Plot regret
        #fig,ax = plt.subplots()
        #ax = plot_regret(res_loaded)
        #plt.tight_layout()
        #plt.savefig(fig_dir + 'regret.pdf')

        print ('Saved skopt figures in ' + fig_dir)
        

def main():

    args=getArgs()
    
    configPath = args.config
    xgb = XGBoostHandler(configPath, args.region)

    if args.skopt_plot: 
        for i in args.fold:
            xgb.skoptPlot(i)
        return

    if args.inputFolder: xgb.setInputFolder(args.inputFolder)
    if args.outputFolder: xgb.setOutputFolder(args.outputFolder)
    if args.params: xgb.setParams(args.params)

    xgb.readData()

    fixRandomSeed(1)

    # looping over the "4 folds"
    for i in args.fold:

        print('===================================================')
        print(f' The #{i} fold of training')
        print('===================================================')

        #xgb.setParams({'eval_metric': ['auc', 'logloss']}, i)
        xgb.set_early_stopping_rounds(20)

        xgb.prepareData(i)

        if args.skopt: 
            try:
                xgb.skoptHP(i)
            except KeyboardInterrupt:
                print('Finishing on SIGINT.')
            continue

        xgb.buildModel(i)

        xgb.trainModel(i)

        xgb.testModel(i)
        print("param: %s, Val AUC: %f" % xgb.getAUC(i))
        print('++++Test AUC++++')
        print(xgb.getAUC(i, sample_set='test') )
        # xgb.plotScore(i, 'test')
        # xgb.plotFeaturesImportance(i)
        xgb.plotROC(i)

        xgb.transformScore(i)

        if args.save: xgb.save(i)

    print('------------------------------------------------------------------------------')
    print('Finished training.')

    return

if __name__ == '__main__':
    main()
