import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--type', action='store', default='weight', help='weight, gain, or cover')
parser.add_argument('-v', '--VBF', action='store_true', help='weight, gain, or cover')
parser.add_argument('-r', '--region', action='store', choices=['zero_jet','one_jet','two_jet','all_jet'], default='all_jet', help='Region to process')
args = parser.parse_args()

import math

#import ROOT
#from root_numpy import root2array

from array import array

import numpy as np


#from sklearn.preprocessing import StandardScaler

#from keras import models
import xgboost as xgb
from xgboost import plot_importance

import matplotlib.pyplot as plt

def process_features():

    region = args.region

    model = xgb.Booster()
    model.load_model('models/%s%s_0.h5' % ('VBF_' if args.VBF else '', region))
    
    xgb.plot_importance(booster =model, importance_type = args.type)


    #plt.rcdefaults()
    fig, ax = plt.subplots()

    importances = model.get_score(importance_type = args.type).values()
    names =  model.get_score(importance_type = args.type).keys()
    print model.get_score(importance_type = args.type).items()

    names = [s.strip('f') for s in names]
    names = [int(i) for i in names]

    print len(importances), " features used to make BDT cuts"

    real_names = [
	     "jet1_pt", "jet1_eta", "jet1_phi*",
             "jet2_pt",	"jet2_eta", "jet2_phi*",
             "met",	"met_phi",
             "dijet_pt", "dijet_eta", "dijet_phi*", "dijet_m",
             "Z_pt", "Z_eta", "dummy", "costheta*"
             ]

    print names
    for i in range(len(names)):
      names[i] = real_names[names[i]]
    print names
    print len(names), " names"

    features = sorted(zip(importances,names))  
    importances = list(reversed([x for x, y in features]))
    names = list(reversed([y for x, y in features]))

    #for i in range(3):
    #    importances.pop()
    #    names.pop()

    y_pos = np.arange(len(names))

    ax.barh(y_pos, importances, align='center',
            color='blue', ecolor='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Importance (%s)' % args.type )
    ax.set_title('Feature Importance Hadronic')

    if not os.path.isdir('features_plots'):
        os.makedirs('features_plots')

    plt.savefig('features_plots/%s_%sfeature.pdf' % (region, 'VBF_' if args.VBF else ''))
    plt.show()


    return


def main():
    process_features()

    return

if __name__ == '__main__':
    main()

