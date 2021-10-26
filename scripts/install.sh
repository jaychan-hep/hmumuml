# create environment from the `environment.yml`
# conda env create -f environment.yml
# conda activate hmumuml

# source /cvmfs/sft.cern.ch/lcg/releases/LCG_97apython3/ROOT/v6.20.06/x86_64-centos7-gcc8-opt/ROOT-env.sh
python3 -m venv hmumumlenv
# include ROOT
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.02/x86_64-centos7-gcc48-opt/bin/thisroot.sh
source hmumumlenv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt


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
