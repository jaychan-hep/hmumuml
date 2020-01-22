#!/usr/bin/env python
import os
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd
from root_pandas import *
import pickle
from tqdm import tqdm
from pdb import set_trace

#pd.options.mode.chained_assignment = None

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', default='data/apply_config.json', help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', default='outputs', help='directory of files containing BDT scores')
    parser.add_argument('-o', '--outputFolder', action='store', default='fitInputs', help='directory for outputs containing category indeces')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'all_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-b', '--nbin', type = int, default = 4, choices = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of BDT bins.')
    parser.add_argument('-v', '--vbf', type = int, default = 0, choices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16], help = 'number of VBF bins.')
    return parser.parse_args()

class ApplyXGBHandler(object):
    "Class for applying XGBoost"

    def __init__(self, configPath, region=''):

        print('===============================')
        print('  ApplyXGBHandler initialized')
        print('===============================')

        self._region = region
        self._inputFolder = ''
        self._inputTree = 'test'
        self._outputFolder = ''
        self._chunksize = 500000

        self.randomIndex = 'eventNumber'

        self.boundaries = [[]]
        self.boundaries_v = [[]]
        self._nbin = 0
        self._nVBF = 0
        self._nggF = 0

        self.readApplyConfig(configPath)

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

        except Exception as e:
            logging.error("Error reading apply configuration '{config}'".format(config=configPath))
            logging.error(e)

    def setBoundaries(self, boundaries, boundaries_v):
        self.boundaries = boundaries
        self.boundaries_v = boundaries_v
        self._nbin = len(boundaries[0]) + len(boundaries_v[0])
        self._nVBF = len(boundaries_v[0])
        self._nggF = len(boundaries[0])

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def applyBDT(self, category):

        outputContainer = self._outputFolder + '/' + self._region
        output_path = outputContainer + '/%s.root' % category
        if not os.path.isdir(outputContainer): os.makedirs(outputContainer)
        if os.path.isfile(output_path): os.remove(output_path)

        f_list = []
        input_path = self._inputFolder + '/' + self._region + '/%s.root' % category

        for data in tqdm(read_root(input_path, key=self._inputTree, chunksize=self._chunksize), desc='XGB INFO: Assigning category index to %s samples' % category, ncols=100):

            out_data = pd.DataFrame()

            for i in range(len(self.boundaries)):

                data_s = data[data[self.randomIndex]%len(self.boundaries) == i]

                if self.boundaries[i][0] != 0:
                    data_s = data_s[data_s['bdt_score_t'] >= self.boundaries[i][0]]
                if self._nVBF == 0:
                    data_s['bdt_category'] = sum([np.heaviside(b-data_s['bdt_score_t'],0) for b in self.boundaries[i]]) + 1
                else:
                    data_s['bdt_category'] = np.heaviside(self.boundaries_v[i][0] - data_s['bdt_score_VBF_t'], 0) * (sum([np.heaviside(b - data_s['bdt_score_t'], 0) for b in self.boundaries[i]])) + sum([np.heaviside(b - data_s['bdt_score_VBF_t'], 0) for b in self.boundaries_v[i]]) + 1

                out_data = pd.concat([out_data, data_s], ignore_index=True, sort=False)

            out_data = out_data.astype({"bdt_category": int})
            out_data.to_root(output_path, key='test', mode='a', index=False)


def main():

    args=getArgs()
   
    nvbf, nggf = args.vbf, args.nbin
    with open('significances/%s/%s%d.json'%(args.region, '%d_'%nvbf if nvbf != 0 else '', nggf)) as f:
        data = json.load(f)

    boundaries_v = data['boundaries_VBF_values'] if nvbf != 0 else [[]]
    boundaries = data['boundaries_values']
 
    configPath = args.config
    xgb = ApplyXGBHandler(configPath, args.region)

    xgb.setInputFolder(args.inputFolder)
    xgb.setOutputFolder(args.outputFolder)

    xgb.setBoundaries(boundaries, boundaries_v)

    with open('data/inputs_config.json') as f:
        config = json.load(f)
    sample_list = config['sample_list']

    for category in sample_list:
        xgb.applyBDT(category)

    os.system('hadd {folder}/{region}/sig.root {folder}/{region}/{{VBF.root,ggF.root,ttH.root,VH.root}}'.format(folder = args.outputFolder, region = args.region))

    return

if __name__ == '__main__':
    main()
