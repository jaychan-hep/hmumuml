#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.13.2018
#
#
#
#
#

import os
from datetime import datetime
from argparse import ArgumentParser
from ttHyy.condor import condor_booklist

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for Hmumu XGBoost analysis.")
    parser.add_argument('-i', '--inputdir', action='store', default='inputs', help='Directory that contains ntuples.')
    return  parser.parse_args()


def arrange_samplelist(channel,category,inputdir):
    samples=[]
    for filename in os.listdir(inputdir):
        if channel not in filename: continue
        if '.root' not in filename: continue
        for section in ['g', '0', '1', '2', '3']:
            #if section != 'g' and category == 'data': continue
            if category in ['ttbar', 'Z', 'stop', 'diboson'] and section != 'g': continue
            for region in ['two_jet','one_jet','zero_jet']:
                samples.append("%s %s %s %s"% (filename,category,section,region))
    return samples

def main():
    args=getArgs()
    date = datetime.now().strftime("%Y-%m-%d-%H-%M")
    inputdir = args.inputdir
    sample_list = {'data':    ['data15','data16','data17','data18'],
                   'ttbar':   ['410472'],
                   'Z':       ['364100','364101','364102','364103','364104','364105','364106','364107','364108','364109','364110','364111','364112','364113',
                               '366300','366301','366303','366304','366306'],
                   'stop':    ['410644','410645'],
                   'diboson': ['363356','363358','364250','364253','364254'],
                   'ttH':     ['344388'],
                   'VH':      ['345103','345104','345105','345098'],
                   'VBF':     ['345106'],
                   'ggF':     ['345097']}

    condor_list = {}

    for category in sample_list:

        condor_list[category] = condor_booklist('SubmitProcessArrays.sh', 'ProcessArrays', category)
        condor_list[category].initialdir_in_arguments()
        condor_list[category].set_JobFlavour('longlunch')

        for channel in sample_list[category]:

            samples = arrange_samplelist(channel,category,inputdir)
            condor_list[category].add_Argument(samples)

        condor_list[category].summary('Basic')
        condor_list[category].submit()

    
if __name__=='__main__':
    main()

