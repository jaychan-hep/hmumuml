# create environment from the `environment.yml`
conda env create -f environment.yml
conda activate hmumuml

# create virtual environment
#conda create -n hmumuml python=3.7 # --clone base
#conda activate hmumuml

# install packages
#conda config --add channels conda-forge
#conda config --add channels anaconda
#conda install root root_numpy
#conda install h5py scikit-learn matplotlib tabulate xgboost pandas root_pandas tqdm

# setting python path
export PATH="`pwd`/scripts:${PATH}"
export PYTHONPATH="`pwd`/scripts:${PYTHONPATH}"
export PATH="`pwd`/hmumuml:$PATH"
export PYTHONPATH="`pwd`/hmumuml:$PYTHONPATH"
