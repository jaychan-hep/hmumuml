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
    parser.add_argument('-f', '--fold', action='store', type=int, nargs='+', choices=[0, 1, 2, 3], default=[0, 1, 2, 3], help='specify the fold for training')
    parser.add_argument('-a', '--algorithm', action='store', choices=["bdt", "NN"], default="NN", help='Which ML algorithm.')
    parser.add_argument('-m', '--mode', action='store', choices=["skopt", "config"], default="skopt", help='To optimize the hypaeparameters with skopt, or do the usual training with config')
    return  parser.parse_args()

def main():

    args=getArgs()

    CONDA_PREFIX = os.getenv("CONDA_PREFIX").replace("/envs/hmumuml", "")

    createScript(f'scripts/submit_hyperparameter_tuning_{args.algorithm}_{args.mode}.sh', f"""#!/bin/bash

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
echo python scripts/train_{args.algorithm}.py -r $region -f $fold {"--skopt" if args.mode == "skopt" else "--save"}
python scripts/train_{args.algorithm}.py -r $region -f $fold {"--skopt" if args.mode == "skopt" else "--save"}
        """)

    condor_list = condor_booklist(f'scripts/submit_hyperparameter_tuning_{args.algorithm}_{args.mode}.sh', 'hyperparameter_tuning')
    condor_list.initialdir_in_arguments()
    condor_list.set_JobFlavour('nextweek')
    condor_list.set_RequestCpus(4)
    condor_list.set_Request_memory("10 GB")
    condor_list.set_Request_disk("1 GB")

    for region in args.region:

        for fold in args.fold:

            condor_list.add_Argument(f"{region} {fold}")

    condor_list.summary('Basic')
    condor_list.submit()

    
if __name__=='__main__':
    main()

