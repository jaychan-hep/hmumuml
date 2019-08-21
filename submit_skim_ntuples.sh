#!/bin/bash

initdir=$1
input=$2
output=$3

cd $initdir
source ./setup_lxplus.sh
echo python skim_ntuples.py -i $input -o $output
python skim_ntuples.py -i $input -o $output
