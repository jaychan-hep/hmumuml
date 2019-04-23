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
if [ $3 == 'g' ]
    then
    runsection=-1
else
    runsection=$3
fi
region=$4
category=$5

cd $initdir
source ./setup_pcuw.sh
echo python process_arrays.py -n $filename -s $runsection -r $region -c $category
if [[ "$category" == *data* ]]; then
python process_arrays.py -n $filename -s $runsection -r $region -c $category -d
else
python process_arrays.py -n $filename -s $runsection -r $region -c $category
fi
