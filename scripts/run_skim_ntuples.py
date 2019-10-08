#!/usr/bin/env python

import os
from argparse import ArgumentParser
from ROOT import TFile
import json

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for Hmumu XGBoost analysis.")
    parser.add_argument('-i', '--inputdir', action='store', default='inputs', help='Directory that contains ntuples.')
    parser.add_argument('-c', '--check', action = 'store_true', default = False, help = 'Check the completeness')
    parser.add_argument('-s', '--skip', action = 'store_true', default = False, help = 'Skip checking the metadata')
    return  parser.parse_args()


def arrange_cmd(channel,category,inputdir, skip):

    cmd = []

    for filename in os.listdir(inputdir):

        if channel not in filename: continue
        if '.root' not in filename: continue
        if os.path.isfile('skimmed_ntuples/%s/%s' % (category, filename)):

            if skip: continue

            print 'Checking %s...' % filename

            f0 = TFile('inputs/%s' % filename)
            t0 = f0.Get('DiMuonNtuple')
            f1 = TFile('skimmed_ntuples/%s/%s' % (category, filename))
            t1 = f1.Get('inclusive')
            t10 = f1.Get('zero_jet')
            t11 = f1.Get('one_jet')
            t12 = f1.Get('two_jet')

            try:
                if (t0.GetEntries("FinalSelection && Muons_Minv_MuMu_Fsr >= 110 && SampleOverlapWeight && EventWeight_MCCleaning_5") == t1.GetEntries()) and (t1.GetEntries() == t10.GetEntries() + t11.GetEntries() + t12.GetEntries()):
                        continue
                else:
                    print "Events after desired selections don't match!! Please check!! Removing the skimmed file..."
                    #os.remove('skimmed_ntuples/%s/%s' % (category, filename))
            except Exception as e:
                print e

        else:
            print 'Missing %s!!' % filename

        if not os.path.isdir('skimmed_ntuples/%s' % category): os.makedirs('skimmed_ntuples/%s' % category)
        cmd.append("python scripts/skim_ntuples.py -i inputs/%s -o skimmed_ntuples/%s/%s" % (filename, category, filename))

    return cmd

def main():

    args=getArgs()
    inputdir = args.inputdir

    with open('data/inputs_config.json') as f:
        config = json.load(f)

    sample_list = config['sample_list']

    cmds = []

    for category in sample_list:

        for channel in sample_list[category]:

            cmd = arrange_cmd(channel,category,inputdir, args.skip)
            cmds.extend(cmd)

    if not args.check:
        for cmd in cmds:
            print cmd
            os.system(cmd)

    if args.check:
        print '='*20
        print 'Unresolved files: ', len(cmds)
        print '='*20
    
if __name__=='__main__':
    main()

