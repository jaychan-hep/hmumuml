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
import pickle
#from sklearn import model_selection
from sklearn.metrics import roc_curve, auc, confusion_matrix
#from sklearn.model_selection import train_test_split
#from keras.optimizers import RMSprop
#from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import xgboost as xgb
from tabulate import tabulate
#from bayes_opt import BayesianOptimization
import matplotlib.pyplot as plt


def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('--cat', '--categorical', action='store_true', help='Create categorical model')
    parser.add_argument('-s', '--signal', action='store', type=signal_multiplier, default='5', help='Number of signal events to use in training as a multiple of the number of background events. Value of 0 means use all signal events.')
    parser.add_argument('-f', '--fold', action='store', type=int, choices=[-1,0,1,2,3], default=-1, help='specify the fold for training')
    parser.add_argument('-n', '--name', action='store', default='test', help='Name of output plot.')
    parser.add_argument('-p', '--params', action='store', type=dict, default=None, help='json string.') #type=json.loads
    parser.add_argument('--numRound', action='store', type=int, default=10000, help='Number of rounds.')
    parser.add_argument('--save', action='store_true', help='Save model weights to HDF5 file')
    parser.add_argument('--sf',type=float,default=1.,help='Scale factor of signal to make the signal and background effectively same in size.')
    parser.add_argument('--roc', action='store_true', help='Plot ROC')
    parser.add_argument('--VBF', action='store_true', help='Train VBF against ggH + data sideband')
    return parser.parse_args()

def signal_multiplier(s):
    s = float(s)
    if s < 0.0:
        raise argparse.ArgumentTypeError('%r must be >= 0!' % s)
    return s

def plotROC(y_train, weight_train, score_train, y_val, weight_val, score_val, y_test, weight_test, score_test, filename, show=False):
    
    fpr_train, tpr_train, _ = roc_curve(y_train, score_train, sample_weight=weight_train)
    roc_auc_train = auc(fpr_train, tpr_train, reorder=True)
    fpr_val, tpr_val, _ = roc_curve(y_val, score_val, sample_weight=weight_val)
    roc_auc_val = auc(fpr_val, tpr_val, reorder=True)
    fpr_test, tpr_test, _ = roc_curve(y_test, score_test, sample_weight=weight_test)
    roc_auc_test = auc(fpr_test, tpr_test, reorder=True)

    print('Training ROC AUC = %f' % roc_auc_train)
    print('Validation ROC AUC = %f' % roc_auc_val)
    print('Test ROC AUC = %f' % roc_auc_test)

    fpr_train = 1.0 - fpr_train
    fpr_val = 1.0 - fpr_val
    fpr_test = 1.0 - fpr_test

    plt.grid(color='gray', linestyle='--', linewidth=1)
    plt.plot(tpr_train, fpr_train, label='Train set, area = %0.6f' % roc_auc_train, color='black', linestyle='dotted')
    plt.plot(tpr_val, fpr_val, label='Val set, area = %0.6f' % roc_auc_val, color='blue', linestyle='dashdot')
    plt.plot(tpr_test, fpr_test, label='Test set, area = %0.6f' % roc_auc_test, color='red', linestyle='dashed')
    plt.plot([0, 1], [1, 0], linestyle='--', color='black', label='Luck')
    plt.xlabel('Signal acceptance')
    plt.ylabel('Background rejection')
    plt.title('Receiver operating characteristic')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xticks(np.arange(0, 1, 0.1))
    plt.yticks(np.arange(0, 1, 0.1))
    plt.legend(loc='lower left', framealpha=1.0)

    #plt.savefig('plots/' + filename + '.png')
    plt.savefig('bdt_plots/' + filename + '.eps')

    if show: plt.show()

    plt.clf()
    return

def getAUC(y_test, weight, score):
    fpr, tpr, _ = roc_curve(y_test, score, sample_weight=weight)
    roc_auc = auc(fpr, tpr, reorder=True)
    return roc_auc

def load_npy(catagories,region):

    badfile=False
    # signal loading
    train_npy=[[],[],[],[]]
    train_wt=[[],[],[],[]]

    for category in catagories:
        print 'Loading %s data...' %(category)
        filedir='arrays/'+category
        for npy in sorted(os.listdir(filedir)):
            if not npy.endswith('_jets.npy'): continue
            if not region in npy: continue
            if '_-1_' in npy: continue
            print npy

            npytemp=[]
            jets=np.load(filedir+'/'+npy)
            if [jets.shape[1],jets.shape[2]]!=[2,5]:
                print 'ERROR: The dimension of jets npy file is not correct!!'
                exit()
            npytemp.append(jets.reshape(len(jets),-1))
            met=np.load(filedir+'/'+npy.replace('jets','met'))
            if [met.shape[1],met.shape[2]]!=[1,2]:
                print 'ERROR: The dimension of met npy file is not correct!!'
                exit()
            npytemp.append(met.reshape(len(met),-1))
            dijet=np.load(filedir+'/'+npy.replace('jets','dijet'))
            if [dijet.shape[1],dijet.shape[2]]!=[1,4]:
                print 'ERROR: The dimension of dijet npy file is not correct!!'
                exit()
            npytemp.append(dijet.reshape(len(dijet),-1))
            dimuon=np.load(filedir+'/'+npy.replace('jets','dimuon'))
            if [dimuon.shape[1],dimuon.shape[2]]!=[1,4]:
                print 'ERROR: The dimension of dimuon npy file is not correct!!'
                exit()
            npytemp.append(dimuon.reshape(len(dimuon),-1))
            npytemp=np.concatenate(npytemp,axis=1)

            wt=np.load(filedir+'/'+npy.replace('jets','weight'))

            train_npy[int(npy.split('_')[4])].append(npytemp)
            train_wt[int(npy.split('_')[4])].append(wt)

    for i in range(4):
        train_npy[i] = np.concatenate(train_npy[i])
        train_wt[i] = np.concatenate(train_wt[i])

    return train_npy, train_wt


def train_model(params, args, region, sigs, bkgs):
    

    print 'Start preparing training samples...'
    # loading data
    train_sig_all, train_sig_wt_all = load_npy(sigs,region)
    train_bkg_all, train_bkg_wt_all = load_npy(bkgs,region)

    SF=args.sf
    num_round=args.numRound

    for i in range(4):

        if args.fold != -1:
            if args.fold != i:
                print '==================================================='
                print 'INFO: skip #%d fold'%i
                print '==================================================='
                continue

        print '==================================================='
        print ' The #%d fold of training' %i
        print '==================================================='

        train_sig = train_sig_all[:]
        train_sig_wt = train_sig_wt_all[:]
        train_bkg = train_bkg_all[:]
        train_bkg_wt = train_bkg_wt_all[:]
        test_sig = train_sig.pop(i)
        test_sig_wt = train_sig_wt.pop(i)
        test_bkg = train_bkg.pop(i)
        test_bkg_wt = train_bkg_wt.pop(i)
        val_sig = train_sig.pop(i%3)
        val_sig_wt = train_sig_wt.pop(i%3)
        val_bkg = train_bkg.pop(i%3)
        val_bkg_wt = train_bkg_wt.pop(i%3)
        train_sig = np.concatenate(train_sig)
        train_sig_wt = np.concatenate(train_sig_wt)
        train_bkg = np.concatenate(train_bkg)
        train_bkg_wt = np.concatenate(train_bkg_wt)

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
        if args.roc: test  = np.concatenate((test_sig  , test_bkg  ))

        y_train_cat = np.concatenate((np.zeros(len(train_sig), dtype=np.uint8), np.ones(len(train_bkg), dtype=np.uint8)))
        y_val_cat   = np.concatenate((np.zeros(len(val_sig)  , dtype=np.uint8), np.ones(len(val_bkg)  , dtype=np.uint8)))
        if args.roc: y_test_cat   = np.concatenate((np.zeros(len(test_sig)  , dtype=np.uint8), np.ones(len(test_bkg)  , dtype=np.uint8)))

        if SF == -1: SF = 1.*len(train_sig_wt)/len(train_bkg_wt)

        y_train_weight = np.concatenate((train_sig_wt*(len(train_sig_wt)+len(train_bkg_wt))*SF/(np.average(train_sig_wt)*(1+SF)*len(train_sig_wt)), train_bkg_wt*(len(train_sig_wt)+len(train_bkg_wt))/(np.average(train_bkg_wt)*(1+SF)*len(train_bkg_wt))))
        y_val_weight   = np.concatenate((val_sig_wt*(len(train_sig_wt)+len(train_bkg_wt))*SF/(np.average(train_sig_wt)*(1+SF)*len(train_sig_wt)), val_bkg_wt*(len(train_sig_wt)+len(train_bkg_wt))/(np.average(train_bkg_wt)*(1+SF)*len(train_bkg_wt))))
        if args.roc: y_test_weight = np.concatenate((test_sig_wt, test_bkg_wt))

        dTrain = xgb.DMatrix(train, label=y_train_cat, weight=y_train_weight)
        dVal   = xgb.DMatrix(val, label=y_val_cat, weight=y_val_weight)
        if args.roc: dTest  = xgb.DMatrix(test, label=y_test_cat, weight=y_test_weight) 
        dTest_sig = xgb.DMatrix(test_sig)

        # train model
        print('Train model.')

        param = {'eval_metric': ['auc'], 'objective': 'binary:logistic'}

        if not args.VBF:
            if args.region == 'two_jet':
                if i == 0: param = {'colsample_bytree': 0.7257488197556066, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 9.761270416292241, 'nthread': 4, 'min_child_weight': 70, 'subsample': 0.8681707406600842, 'eta': 0.01587391341501056, 'max_bin': 174, 'objective': 'binary:logistic', 'alpha': 0.5757150887260754, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 9.55668109157072, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.8722973038051791, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.09021877137744, 'nthread': 4, 'min_child_weight': 91, 'subsample': 0.9054450896669136, 'eta': 0.011252584084082026, 'max_bin': 498, 'objective': 'binary:logistic', 'alpha': 0.7689816483091817, 'tree_method': 'hist', 'max_depth': 11, 'gamma': 8.684593658170515, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.8925064691181133, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 11.579624742406015, 'nthread': 4, 'min_child_weight': 76, 'subsample': 0.864807274924633, 'eta': 0.013262139226209148, 'max_bin': 433, 'objective': 'binary:logistic', 'alpha': 0.979241365658599, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 2.927135805482568, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.8167904025197066, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.793368112649453, 'nthread': 4, 'min_child_weight': 99, 'subsample': 0.8383593536636893, 'eta': 0.016262958003538067, 'max_bin': 498, 'objective': 'binary:logistic', 'alpha': 0.9076947438308036, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 9.8404132158743, 'booster': 'gbtree'}
            if args.region == 'one_jet':
                if i == 0: param = {'colsample_bytree': 0.8516875525290168, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.073469726700036, 'nthread': 4, 'min_child_weight': 95, 'subsample': 0.8835434038296405, 'eta': 0.04567924091321862, 'max_bin': 159, 'objective': 'binary:logistic', 'alpha': 0.5853037734137767, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 9.34012026138051, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.7473310500984237, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.220424562681167, 'nthread': 4, 'min_child_weight': 72, 'subsample': 0.8460636683678937, 'eta': 0.05709245411458616, 'max_bin': 209, 'objective': 'binary:logistic', 'alpha': 0.5655954125660099, 'tree_method': 'hist', 'max_depth': 7, 'gamma': 9.962826948566118, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.7476582376917866, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 5.072576391911285, 'nthread': 4, 'min_child_weight': 83, 'subsample': 0.9255645266985586, 'eta': 0.012655139250546149, 'max_bin': 175, 'objective': 'binary:logistic', 'alpha': 0.9817626015700198, 'tree_method': 'hist', 'max_depth': 8, 'gamma': 9.698398119720297, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.9411793751304067, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.339041745329173, 'nthread': 4, 'min_child_weight': 72, 'subsample': 0.8416196485456633, 'eta': 0.036963489541875734, 'max_bin': 483, 'objective': 'binary:logistic', 'alpha': 0.7318121711587071, 'tree_method': 'hist', 'max_depth': 6, 'gamma': 9.920877239156704, 'booster': 'gbtree'}
            if args.region == 'zero_jet':
                if i == 0: param = {'colsample_bytree': 0.8426641441504055, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 17.85280625044496, 'nthread': 4, 'min_child_weight': 8, 'subsample': 0.8760923929463384, 'eta': 0.02425678476573181, 'max_bin': 493, 'objective': 'binary:logistic', 'alpha': 0.9203532274195957, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 8.596262929835477, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.9421110595840385, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 6.613581929564984, 'nthread': 4, 'min_child_weight': 0, 'subsample': 0.888675974241132, 'eta': 0.08834450850226566, 'max_bin': 479, 'objective': 'binary:logistic', 'alpha': 0.6870495525543625, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 8.941953692455638, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.7468393982398718, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.94548193205644, 'nthread': 4, 'min_child_weight': 0, 'subsample': 0.8286254163297746, 'eta': 0.025996597138504364, 'max_bin': 497, 'objective': 'binary:logistic', 'alpha': 0.8075286032632472, 'tree_method': 'hist', 'max_depth': 6, 'gamma': 5.286967312650763, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.9048460274211345, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 18.977173059816145, 'nthread': 4, 'min_child_weight': 3, 'subsample': 0.950624734630124, 'eta': 0.08833744590953348, 'max_bin': 493, 'objective': 'binary:logistic', 'alpha': 0.8706841149388309, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 0.9305799494143896, 'booster': 'gbtree'}

        else:
            if args.region == 'two_jet':
                if i == 0: param = {'colsample_bytree': 0.7043186973028924, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 17.534207040049804, 'nthread': 4, 'min_child_weight': 84, 'subsample': 0.9015488912547408, 'eta': 0.013022638140091152, 'max_bin': 388, 'objective': 'binary:logistic', 'alpha': 0.6458425536142127, 'tree_method': 'hist', 'max_depth': 16, 'gamma': 9.05740060024965, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.7648457680143608, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.642865254752067, 'nthread': 4, 'min_child_weight': 84, 'subsample': 0.9573743997775728, 'eta': 0.025621698652061584, 'max_bin': 242, 'objective': 'binary:logistic', 'alpha': 0.6302303501990663, 'tree_method': 'hist', 'max_depth': 17, 'gamma': 6.874062235810102, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.7038591928800596, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 7.190725611435425, 'nthread': 4, 'min_child_weight': 87, 'subsample': 0.857531357296296, 'eta': 0.010455729603020603, 'max_bin': 190, 'objective': 'binary:logistic', 'alpha': 0.5453194832154675, 'tree_method': 'hist', 'max_depth': 18, 'gamma': 8.842040464231504, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.7291461442166165, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 17.408239007167698, 'nthread': 4, 'min_child_weight': 86, 'subsample': 0.8338040537780992, 'eta': 0.011940007263631038, 'max_bin': 97, 'objective': 'binary:logistic', 'alpha': 0.7360031431363234, 'tree_method': 'hist', 'max_depth': 21, 'gamma': 7.026494649811172, 'booster': 'gbtree'}

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
        if args.roc:
            score_test  = bst.predict(dTest)
            score_train = bst.predict(dTrain)
        score_test_sig  = bst.predict(dTest_sig)

        #scoresig=bst.predict(xgb.DMatrix(val_sig))
        #np.save('score.npy',scoresig)
        aucValue = getAUC(y_val_cat, y_val_weight, score)
        print("param: %s, Val AUC: %s" % (param, aucValue))

        if args.roc:
            plotROC(y_train_cat, y_train_weight, score_train, y_val_cat, y_val_weight, score, y_test_cat, y_test_weight, score_test, 'roc_%s_%d'%(region, i), show=True)

        tsf = QuantileTransformer(n_quantiles=1000, output_distribution='uniform', subsample=1000000000, random_state=0)
        #plt.hist(score_test_sig, bins='auto')
        #plt.show()
        tsf.fit(score_test_sig.reshape(-1, 1))
        score_test_sig_t=tsf.transform(score_test_sig.reshape(-1, 1)).reshape(-1)
        #plt.hist(score_test_sig_t, bins='auto')
        #plt.show()

        # save model
        if not os.path.isdir('models'):
            os.makedirs('models')
        if args.save:
            print('Saving model')
            bst.save_model('models/%s%s_%d.h5' % ('' if not args.VBF else 'VBF_',region,i))

            with open('models/score_transformer_%s%s_%d.pkl'%(region, '_VBF' if args.VBF else '', i), 'wb') as f:
                pickle.dump(tsf, f, -1)

    return True




def main():
    args=getArgs()
    params=args.params
    

    print '=============================================================================='
    sigs=['VBF', 'ggF'] if not args.VBF else ['VBF']
    
    print 'INFO:  Training as signal on:  ', sigs

    bkgs = ['data'] if not args.VBF else ['data']
    print 'INFO:  Traing as backgroun on:  ', bkgs
    print '------------------------------------------------------------------------------'

    while not train_model(params, args, args.region, sigs, bkgs):
        os.system('python RecoverMissingFiles.py -b')

    print '------------------------------------------------------------------------------'
    print 'Finished training.'

    return

if __name__ == '__main__':
    main()
