#!/bin/bash
#
#
#
#    Created by Jay Chan
#
#        10.13.2018
#
#
#
#
initdir=$1
filename=$2
category=$3
if [[ $4 == "g" ]]; then
    runsection=-1
else
    runsection=$4
fi
region=$5


cd $initdir
source ./setup_lxplus.sh
echo python process_arrays.py -n $filename -s $runsection -r $region -c $category
if [[ $category == *"data"* ]]; then
    echo "data!"
    python process_arrays.py -n $filename -s $runsection -r $region -c $category -d
else
    echo "MC!"
    python process_arrays.py -n $filename -s $runsection -r $region -c $category
fi
