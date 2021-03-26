#!/usr/bin/env python

import os
from datetime import datetime
from argparse import ArgumentParser
from condor import condor_booklist, createScript
import json

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser(description="Hyperparameter tuning.")
    parser.add_argument('-r', '--region', action='store', nargs="+", default=["zero_jet", "one_jet", "two_jet", "VBF"], help='Regions to run.')
    return  parser.parse_args()

def main():

    args=getArgs()

    CONDA_PREFIX = os.getenv("CONDA_PREFIX").replace("/envs/hmumuml", "")

    createScript('scripts/submit_hyperparameter_tuning.sh', f"""#!/bin/bash

initdir=$1
region=$2
fold=$3

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('{CONDA_PREFIX}/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "{CONDA_PREFIX}/etc/profile.d/conda.sh" ]; then
        . "{CONDA_PREFIX}/etc/profile.d/conda.sh"
    else
        export PATH="{CONDA_PREFIX}/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

cd $initdir
source scripts/setup.sh
echo python scripts/train_bdt.py -r $region -f $fold --skopt
python scripts/train_bdt.py -r $region -f $fold --skopt
        """)

    condor_list = condor_booklist('scripts/submit_hyperparameter_tuning.sh', 'hyperparameter_tuning')
    condor_list.initialdir_in_arguments()
    condor_list.set_JobFlavour('testmatch')
    condor_list.set_RequestCpus(4)
    condor_list.set_Request_memory("10 GB")
    condor_list.set_Request_disk("1 GB")

    for region in args.region:

        for fold in ["0", "1", "2", "3"]:

            condor_list.add_Argument(f"{region} {fold}")

    condor_list.summary('Basic')
    condor_list.submit()

    
if __name__=='__main__':
    main()

