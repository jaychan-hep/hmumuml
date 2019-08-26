#!/bin/bash
#
#
#
#    Created by Jay Chan
#
#        3.10.2019
#
#
#
#
initdir=$1
vb=$2
nbin=$3
nvbf=$4


cd $initdir
source ./setup_lxplus.sh
echo python scripts/categorization_optimization_2D.py -t --vb $vb
python scripts/categorization_optimization_2D.py -t --vb $vb -b $nbin -v $nvbf -r two_jet -f 1 #--floatB
