#!/bin/bash
# include ROOT
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.02/x86_64-centos7-gcc48-opt/bin/thisroot.sh
source hmumumlenv/bin/activate
#export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
#source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
# conda activate hmumumlenv

export PATH="`pwd`:${PATH}"
export PYTHONPATH="`pwd`:${PYTHONPATH}"
export THEANO_FLAGS="gcc.cxxflags='-march=core2'"

export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"

export PATH="`pwd`/hmumuml:$PATH"
export PYTHONPATH="`pwd`/hmumuml:$PYTHONPATH"
