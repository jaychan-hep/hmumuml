#!/bin/bash
#
#
#
#    Created by Jay Chan
#
#        8.13.2018
#
#
#
#
initdir=$1
filename=$2
category=$3


cd $initdir
source ./setup_lxplus.sh
echo python process_arrays.py -n $filename $runsection -r $region -c $category
if [[ $category == *"data"* ]]; then
    echo "data!"
    python process_arrays.py -n $filename -t -r two_jet -c $category -d
    python process_arrays.py -n $filename -v -r two_jet -c $category -d
    python process_arrays.py -n $filename -r two_jet -c $category -d
    python process_arrays.py -n $filename -t -r one_jet -c $category -d
    python process_arrays.py -n $filename -v -r one_jet -c $category -d
    python process_arrays.py -n $filename -r one_jet -c $category -d
    python process_arrays.py -n $filename -t -r zero_jet -c $category -d
    python process_arrays.py -n $filename -v -r zero_jet -c $category -d
    python process_arrays.py -n $filename -r zero_jet -c $category -d
else
    echo "MC!"
    python process_arrays.py -n $filename -t -r two_jet -c $category
    python process_arrays.py -n $filename -v -r two_jet -c $category
    python process_arrays.py -n $filename -r two_jet -c $category
    python process_arrays.py -n $filename -t -r one_jet -c $category
    python process_arrays.py -n $filename -v -r one_jet -c $category
    python process_arrays.py -n $filename -r one_jet -c $category
    python process_arrays.py -n $filename -t -r zero_jet -c $category
    python process_arrays.py -n $filename -v -r zero_jet -c $category
    python process_arrays.py -n $filename -r zero_jet -c $category
fi
