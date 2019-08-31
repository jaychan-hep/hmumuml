#!/bin/bash

initdir=$1
input=$2
output=$3

cd $initdir
source scripts/setup.sh
echo python scripts/skim_ntuples.py -i $input -o $output
python scripts/skim_ntuples.py -i $input -o $output
