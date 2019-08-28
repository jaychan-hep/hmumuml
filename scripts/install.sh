# setup environment
setupATLAS
lsetup "root  6.12.06-x86_64-slc6-gcc62-opt"
lsetup python
virtualenv --python=python2.7 ve
source ve/bin/activate

# install packages
pip install pip --upgrade
pip install h5py sklearn matplotlib tabulate xgboost pandas root_pandas progressbar
pip install --upgrade https://github.com/rootpy/root_numpy/zipball/master

# setting python path
export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"
export PATH="`pwd`/hmumuml:$PATH"
export PYTHONPATH="`pwd`/hmumuml:$PYTHONPATH"
