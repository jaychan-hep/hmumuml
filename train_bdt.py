#!/usr/bin/env python
import os
from argparse import ArgumentParser
#import json
#import traceback
#import math
import ROOT
#from root_numpy import root2array
import numpy as np
from ttHyy import *
#import pickle
#from sklearn import model_selection
from sklearn.metrics import roc_curve, auc, confusion_matrix
#from sklearn.model_selection import train_test_split
#from keras.optimizers import RMSprop
#from keras.callbacks import EarlyStopping, ModelCheckpoint
#from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from tabulate import tabulate
#from bayes_opt import BayesianOptimization


def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('--cat', '--categorical', action='store_true', help='Create categorical model')
    parser.add_argument('-s', '--signal', action='store', type=signal_multiplier, default='5', help='Number of signal events to use in training as a multiple of the number of background events. Value of 0 means use all signal events.')
    parser.add_argument('-n', '--name', action='store', default='test', help='Name of output plot.')
    parser.add_argument('-p', '--params', action='store', type=dict, default=None, help='json string.') #type=json.loads
    parser.add_argument('--numRound', action='store', type=int, default=10000, help='Number of rounds.')
    parser.add_argument('--save', action='store_true', help='Save model weights to HDF5 file')
    parser.add_argument('--sf',type=float,default=1.,help='Scale factor of signal to make the signal and background effectively same in size.')
    return parser.parse_args()

def signal_multiplier(s):
    s = float(s)
    if s < 0.0:
        raise argparse.ArgumentTypeError('%r must be >= 0!' % s)
    return s


def getAUC(y_test, score):
    fpr, tpr, _ = roc_curve(y_test, score)
    roc_auc = auc(fpr, tpr)
    return roc_auc

def load_npy(catagories,region):

    badfile=False
    # signal loading
    train_npy=[]
    val_npy=[]
    train_wt=[]
    val_wt=[]
    for category in catagories:
        print 'Loading %s data...' %(category)
        filedir='arrays/'+category
        for npy in sorted(os.listdir(filedir)):
            if not npy.endswith('_%s_jets.npy' %(region)):
                continue
            if not ('_val_' in npy or '_train_' in npy):
                continue
            print npy
            try:
                npytemp=[]
                jets=np.load(filedir+'/'+npy)
                if [jets.shape[1],jets.shape[2]]!=[2,5]:
                    print 'ERROR: The dimension of jets npy file is not correct!!' 
                    error
                npytemp.append(jets.reshape(len(jets),-1))
                met=np.load(filedir+'/'+npy.replace('jets','met'))
                if [met.shape[1],met.shape[2]]!=[1,2]:
                    print 'ERROR: The dimension of met npy file is not correct!!'
                    error
                npytemp.append(met.reshape(len(met),-1))
                muons=np.load(filedir+'/'+npy.replace('jets','muons'))
                if [muons.shape[1],muons.shape[2]]!=[2,4]:
                    print 'ERROR: The dimension of muons npy file is not correct!!'
                    error
                npytemp.append(muons.reshape(len(muons),-1))
                dijet=np.load(filedir+'/'+npy.replace('jets','dijet'))
                if [dijet.shape[1],dijet.shape[2]]!=[1,6]:
                    print 'ERROR: The dimension of dijet npy file is not correct!!'
                    error
                npytemp.append(dijet.reshape(len(dijet),-1))
                dimuon=np.load(filedir+'/'+npy.replace('jets','dimuon'))
                if [dimuon.shape[1],dimuon.shape[2]]!=[1,5]:
                    print 'ERROR: The dimension of dimuon npy file is not correct!!'
                    error
                npytemp.append(dimuon.reshape(len(dimuon),-1))
                extras=np.load(filedir+'/'+npy.replace('jets','extras'))
                if [extras.shape[1],extras.shape[2]]!=[1,6]:
                    print 'ERROR: The dimension of extras npy file is not correct!!'
                    error
                npytemp.append(extras.reshape(len(extras),-1))
                npytemp=np.concatenate(npytemp,axis=1)
                #npytemp=utils.restrictSample(npytemp, len(npytemp), 1)

                wt=np.load(filedir+'/'+npy.replace('jets','weight'))

                if '_train_' in npy:
                    train_npy.append(npytemp)
                    train_wt.append(wt)
                else:
                    val_npy.append(npytemp)
                    val_wt.append(wt)
            except:
                badfile=True
                print('ERROR: Unable to load \"'+ filedir+'/'+npy+'\"')
                for filename in os.listdir('inputs'):
                    if npy.split('_')[1] in filename:
                        print 'Original file: \"%s\"'%(filename)
                        try:
                            if txt: pass
                        except:
                            txt = open('arrays/badsamplelist.txt','w')
                        if '_train_' in npy: section='_train'
                        elif '_val_' in npy: section='_val'
                        else: section=''
                        txt.write('%s %s %s %s badsample\n'%(filename, section, region, category))
                        break
                return 0, 0, 0, 0, badfile
    train=np.concatenate(train_npy)
    val=np.concatenate(val_npy)
    train_wt=np.concatenate(train_wt)
    val_wt=np.concatenate(val_wt)

    return train, val, train_wt, val_wt, badfile

def train_model(params, args, region, sigs, bkgs):
    

    print 'Start preparing training samples...'
    # loading data
    train_sig, val_sig, train_sig_wt, val_sig_wt, badfilesig = load_npy(sigs,region)
    train_bkg, val_bkg, train_bkg_wt, val_bkg_wt, badfilebkg = load_npy(bkgs,region)

    #train_sig = utils.restrictSample(train_sig, len(train_bkg), args.signal)
    #val_sig = utils.restrictSample(val_sig, len(val_bkg), args.signal)
    #train_sig_wt = utils.restrictSample(train_sig_wt, len(train_bkg_wt), args.signal)
    #val_sig_wt = utils.restrictSample(val_sig_wt, len(val_bkg_wt), args.signal)

    if badfilesig or badfilebkg:
        return False

    # split data into train, val, and test samples
    print('Splitting data.')

    headers = ['Sample', 'Total', 'Training', 'Validation']
    sample_size_table = [
        ['Signal'    , len(train_sig)+len(val_sig), len(train_sig), len(val_sig)],
        ['Background', len(train_bkg)+len(val_bkg), len(train_bkg), len(val_bkg)],
    ]
    print tabulate(sample_size_table, headers=headers, tablefmt='simple')

    # organize data for training
    print('Organizing data for training.')

    train = np.concatenate((train_sig, train_bkg))
    val   = np.concatenate((val_sig  , val_bkg  ))

    y_train_cat = np.concatenate((np.zeros(len(train_sig), dtype=np.uint8), np.ones(len(train_bkg), dtype=np.uint8)))
    y_val_cat   = np.concatenate((np.zeros(len(val_sig)  , dtype=np.uint8), np.ones(len(val_bkg)  , dtype=np.uint8)))

    SF=args.sf
    if SF == -1: SF = 1.*len(train_sig_wt)/len(train_bkg_wt)

    y_train_weight = np.concatenate((train_sig_wt*(len(train_sig_wt)+len(train_bkg_wt))*SF/(np.average(train_sig_wt)*(1+SF)*len(train_sig_wt)), train_bkg_wt*(len(train_sig_wt)+len(train_bkg_wt))/(np.average(train_bkg_wt)*(1+SF)*len(train_bkg_wt))))
    y_val_weight   = np.concatenate((val_sig_wt*(len(train_sig_wt)+len(train_bkg_wt))*SF/(np.average(train_sig_wt)*(1+SF)*len(train_sig_wt)), val_bkg_wt*(len(train_sig_wt)+len(train_bkg_wt))/(np.average(train_bkg_wt)*(1+SF)*len(train_bkg_wt))))

    dTrain = xgb.DMatrix(train, label=y_train_cat, weight=y_train_weight)
    dVal = xgb.DMatrix(val, label=y_val_cat, weight=y_val_weight)
    # train model
    print('Train model.')
    num_round=args.numRound

    param = {'eval_metric': ['auc'], 'objective': 'binary:logistic'}
    if args.region == 'two_jet': param = {'colsample_bytree': 0.966403297821083, 'silent': 1, 'eval_metric': ['rmse', 'auc', 'logloss'], 'grow_policy': 'lossguide', 'max_delta_step': 0.2514774307847155, 'min_child_weight': 3.9459278013718415, 'subsample': 0.9429725422180117, 'num_boost_round': 100000, 'eta': 0.10337996243247606, 'max_bin': 364, 'objective': 'binary:logistic', 'alpha': 0.8467431405313621, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 2.216870614863294, 'booster': 'gbtree'}
    if args.region == 'one_jet': param = {'colsample_bytree': 0.9688574977915417, 'silent': 1, 'eval_metric': ['rmse', 'auc', 'logloss'], 'grow_policy': 'lossguide', 'max_delta_step': 27.31042949483293, 'min_child_weight': 2.6725730039632007, 'subsample': 0.7417275532809493, 'num_boost_round': 100000, 'eta': 0.05172290311989361, 'max_bin': 322, 'objective': 'binary:logistic', 'alpha': 0.3959167077938971, 'tree_method': 'hist', 'max_depth': 8, 'gamma': 0.01937210341151405, 'booster': 'gbtree'}
    if args.region == 'zero_jet': param = {'colsample_bytree': 0.7502402719138584, 'silent': 1, 'eval_metric': ['rmse', 'auc', 'logloss'], 'grow_policy': 'lossguide', 'max_delta_step': 20.5038891025236, 'min_child_weight': 0.22707394602410563, 'subsample': 0.864312471949317, 'num_boost_round': 100000, 'eta': 0.015811643343602133, 'max_bin': 353, 'objective': 'binary:logistic', 'alpha': 0.43668191400379675, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 1.5094587115179945, 'booster': 'gbtree'}

    if params:
        for key in params:
            param[key] = params[key]

    print(param)

    evallist  = [(dTrain, 'train'), (dVal, 'eval')]
    evals_result = {}
    eval_result_history = []
    try:
        bst = xgb.train(param, dTrain, num_round, evals=evallist, early_stopping_rounds=10, evals_result=evals_result)
    except KeyboardInterrupt:
        print('Finishing on SIGINT.')

    # test model
    print('Test model.')

    # will do: split data to have train, val and test datasets
    score = bst.predict(dVal)
    #scoresig=bst.predict(xgb.DMatrix(val_sig))
    #np.save('score.npy',scoresig)
    aucValue = getAUC(y_val_cat, score)
    print("param: %s, Val AUC: %s" % (param, aucValue))

    print(" ")
    print("======reversing training and validation samples======")
    evallist  = [(dVal, 'train'), (dTrain, 'eval')]
    evals_result = {}
    eval_result_history = []
    try:
        bst_r = xgb.train(param, dVal, num_round, evals=evallist, early_stopping_rounds=10, evals_result=evals_result)
    except KeyboardInterrupt:
        print('Finishing on SIGINT.')

    # save model
    if not os.path.isdir('models'):
        os.makedirs('models')
    if args.save:
        print('Saving model')
        signame=''
        #for sig in sigs:
        #    signame+=sig+'_'
        bst.save_model('models/%s%s.h5' % (signame,region))
        bst_r.save_model('models/%s%s_r.h5' % (signame,region))

    return True




def main():
    args=getArgs()
    params=args.params
    

    print '=============================================================================='
    sigs=['VBF', 'ggF']
    
    print 'INFO:  Training as signal on:  ', sigs

    bkgs = ['data']
    print 'INFO:  Traing as backgroun on:  ', bkgs
    print '------------------------------------------------------------------------------'

    if args.region == 'two_jet':
        while not train_model(params, args, 'two_jet',sigs,bkgs):
            os.system('python RecoverMissingFiles.py -b')
    elif args.region == 'one_jet':
        while not train_model(params, args, 'one_jet',sigs,bkgs):
            os.system('python RecoverMissingFiles.py -b')
    else:
        while not train_model(params, args, 'zero_jet',sigs,bkgs):
            os.system('python RecoverMissingFiles.py -b')

    print '------------------------------------------------------------------------------'
    print 'Finished training.'

    return

if __name__ == '__main__':
    main()
