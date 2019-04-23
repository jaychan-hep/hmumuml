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
#source ./setup_lxplus.sh
source ./setup_pcuw.sh
echo python categorization_optimization_2D.py -t --vb $vb
python categorization_optimization_2D.py -t --vb $vb -b $nbin -v $nvbf -r two_jet --floatB
