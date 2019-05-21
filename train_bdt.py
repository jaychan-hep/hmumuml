#!/usr/bin/env python
import os
from argparse import ArgumentParser
#import json
import numpy as np
import pickle
from sklearn.metrics import roc_curve, auc, confusion_matrix
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
                if i == 0: param = {'colsample_bytree': 0.8527012148953288, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 10.553214676316191, 'nthread': 4, 'min_child_weight': 80, 'subsample': 0.828340283247559, 'eta': 0.010370796056718267, 'max_bin': 298, 'objective': 'binary:logistic', 'alpha': 0.7657049978855025, 'tree_method': 'hist', 'max_depth': 12, 'gamma': 8.115786988227082, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.9321915820256219, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 5.360530732225799, 'nthread': 4, 'min_child_weight': 81, 'subsample': 0.8253978830871902, 'eta': 0.010033110917550287, 'max_bin': 484, 'objective': 'binary:logistic', 'alpha': 0.9978578033575904, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 6.742010099411358, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.861418273525663, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 10.225321789206147, 'nthread': 4, 'min_child_weight': 85, 'subsample': 0.8327194425882263, 'eta': 0.01295663946341382, 'max_bin': 463, 'objective': 'binary:logistic', 'alpha': 0.8161676777259822, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 8.641744204860625, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.7377243934838316, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 17.027128054495513, 'nthread': 4, 'min_child_weight': 65, 'subsample': 0.8624032913095777, 'eta': 0.013933978214256798, 'max_bin': 278, 'objective': 'binary:logistic', 'alpha': 0.9678653442793248, 'tree_method': 'hist', 'max_depth': 10, 'gamma': 5.539640288333616, 'booster': 'gbtree'}
            if args.region == 'one_jet':
                if i == 0: param = {'colsample_bytree': 0.7254614242791069, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.76323656339139, 'nthread': 4, 'min_child_weight': 52, 'subsample': 0.8971345545609812, 'eta': 0.024825047692767308, 'max_bin': 498, 'objective': 'binary:logistic', 'alpha': 0.783664733457266, 'tree_method': 'hist', 'max_depth': 8, 'gamma': 9.499025174985997, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.8678248992848352, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 6.379303638241522, 'nthread': 4, 'min_child_weight': 99, 'subsample': 0.8039105350791207, 'eta': 0.02085088986153394, 'max_bin': 417, 'objective': 'binary:logistic', 'alpha': 0.8619028281435293, 'tree_method': 'hist', 'max_depth': 6, 'gamma': 5.342910682916136, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.8992561116076159, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.36288803771862, 'nthread': 4, 'min_child_weight': 47, 'subsample': 0.8076291939603779, 'eta': 0.024416723024429425, 'max_bin': 247, 'objective': 'binary:logistic', 'alpha': 0.5988466651618689, 'tree_method': 'hist', 'max_depth': 6, 'gamma': 9.929192198863149, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.7127953843750405, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 5.012089503301383, 'nthread': 4, 'min_child_weight': 3, 'subsample': 0.8605905092917305, 'eta': 0.05357610617513063, 'max_bin': 214, 'objective': 'binary:logistic', 'alpha': 0.6048464152433651, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 9.347408526019244, 'booster': 'gbtree'}
            if args.region == 'zero_jet':
                if i == 0: param = {'colsample_bytree': 0.683271075285466, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 15.321525796631692, 'nthread': 4, 'min_child_weight': 2, 'subsample': 0.9450044173045702, 'eta': 0.050669687647475536, 'max_bin': 463, 'objective': 'binary:logistic', 'alpha': 0.7913567652445119, 'tree_method': 'hist', 'max_depth': 6, 'gamma': 9.438527347376906, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.5561480915131072, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 5.848464619299429, 'nthread': 4, 'min_child_weight': 6, 'subsample': 0.9349358252932701, 'eta': 0.07970608139168346, 'max_bin': 488, 'objective': 'binary:logistic', 'alpha': 0.8794473874970213, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 9.74750314580846, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.8750321391999238, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.367039456480164, 'nthread': 4, 'min_child_weight': 5, 'subsample': 0.809682361123101, 'eta': 0.04550901517459201, 'max_bin': 497, 'objective': 'binary:logistic', 'alpha': 0.7381933733395216, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 7.786106495754567, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.9865832250991617, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.12287956956653, 'nthread': 4, 'min_child_weight': 1, 'subsample': 0.823050380776366, 'eta': 0.01182389016000303, 'max_bin': 486, 'objective': 'binary:logistic', 'alpha': 0.5454817416782773, 'tree_method': 'hist', 'max_depth': 5, 'gamma': 0.08280285839535842, 'booster': 'gbtree'}

        else:
            if args.region == 'two_jet':
                if i == 0: param = {'colsample_bytree': 0.7078849377895996, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 18.171924379375735, 'nthread': 4, 'min_child_weight': 94, 'subsample': 0.8691980661542644, 'eta': 0.02318341345523771, 'max_bin': 132, 'objective': 'binary:logistic', 'alpha': 0.6773602112327002, 'tree_method': 'hist', 'max_depth': 14, 'gamma': 7.6334636286505635, 'booster': 'gbtree'}
                if i == 1: param = {'colsample_bytree': 0.8822897622218241, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 19.723428301439213, 'nthread': 4, 'min_child_weight': 51, 'subsample': 0.8752836862333582, 'eta': 0.011483290580592843, 'max_bin': 85, 'objective': 'binary:logistic', 'alpha': 0.8631670426284296, 'tree_method': 'hist', 'max_depth': 12, 'gamma': 8.304498702340316, 'booster': 'gbtree'}
                if i == 2: param = {'colsample_bytree': 0.7680078853963639, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 6.423878577794825, 'nthread': 4, 'min_child_weight': 64, 'subsample': 0.8514580979287965, 'eta': 0.012685274215291155, 'max_bin': 485, 'objective': 'binary:logistic', 'alpha': 0.7171962736450721, 'tree_method': 'hist', 'max_depth': 14, 'gamma': 9.801379495128236, 'booster': 'gbtree'}
                if i == 3: param = {'colsample_bytree': 0.7357819073064397, 'silent': 1, 'eval_metric': ['logloss', 'auc'], 'grow_policy': 'lossguide', 'max_delta_step': 12.171231108096954, 'nthread': 4, 'min_child_weight': 61, 'subsample': 0.8525236587842495, 'eta': 0.011323242640265092, 'max_bin': 285, 'objective': 'binary:logistic', 'alpha': 0.9818150563860808, 'tree_method': 'hist', 'max_depth': 22, 'gamma': 8.693854178126491, 'booster': 'gbtree'}

        if params:
            for key in params:
                param[key] = params[key]

        #param = {}
        param['eval_metric'] = ['auc', 'logloss']
        #param['eval_metric'] = ['logloss', 'auc']
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

        # transform the scores
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
