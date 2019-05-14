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
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for MonoHbb XGBoost analysis.")
    parser.add_argument('-i', '--inputdir', action='store', default='inputs', help='Directory that contains ntuples.')
    return  parser.parse_args()


def main():
    args=getArgs()
    date = datetime.now().strftime("%Y-%m-%d-%H-%M")
    inputdir=args.inputdir
    sample_list=['data','ttH','VH','VBF','ggF','Z','ttbar','stop','diboson']

    condor_list={}
    for channel in sample_list:
        #if channel=='Znunu': continue
        condor_list[channel] = condor_booklist('SubmitRecoverMissingFiles.sh', 'RecoverProcessArrays', channel)
        condor_list[channel].initialdir_in_arguments()
        condor_list[channel].set_JobFlavour('longlunch')

    text_file = open("arrays/missing_samplelist.txt", "r")
    lines = text_file.read().splitlines()
    for line in lines:
        args=str(line).split(' ')
        if args[1] == '-1':
            section='g'
        else:
            section= args[1]
        
        #print line
        condor_list[args[3]].add_Argument("%s %s %s %s"% (args[0],section,args[2],args[3]))
        #break
        
    for channel in condor_list:
        condor_list[channel].submit()


    
if __name__=='__main__':
    main()
