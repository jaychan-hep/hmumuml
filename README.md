# Welcome to HmumuML

HmumuML is a package that can be used to study the XGBoost BDT categorization in Hmumu analysis.

## Setup

### First time setup on lxplus

First checkout the code:

```
git clone ssh://git@gitlab.cern.ch:7999/wisc_atlas/hmumuml.git [-b your_branch]
cd HmumuML
```

Then, setup the virtualenv:

```
setupATLAS
lsetup "root  6.12.06-x86_64-slc6-gcc62-opt"
lsetup python
virtualenv --python=python2.7 ve
source ve/bin/activate
```

Then, checkout necessary packages:

```
pip install pip --upgrade
pip install h5py sklearn matplotlib tabulate xgboost
pip install --upgrade https://github.com/rootpy/root_numpy/zipball/master
```

### Normal setup on lxplus

After setting up the environment for the first time, you can return to this setup by doing `source setup_lxplus.sh`

## Scripts to run the tasks

### Producing training inputs and reduced root-files which the trained model will be applied to

The core script is `process_arrays.py`. The script will apply the given preselection cuts to input samples and produce corresponding ntuples, numpy arrays files, and a text file containing the information of averaged weight of the sample. You can directly run it for any specified input files. Since there are too many samples in H->mumu, it is more convinient and time-efficient to use `SubmitProcessArrays.py` which will find all of the input files and run the core script by submitting the condor jobs.


```
python SubmitProcessArrays.py
```
or
```
python process_arrays.py [-n INPUT_FILE_NAME] [-r REGION] [-s SECTION] [-c CATEGORY] [-d]

Usage:
  -n, --name          The name of the input root-file. The mother folder ('inputs/') shall not be included in the name. 
  -r, --region        The region of interest. Choices: 'zero_jet', 'one_jet', or 'two_jet'.
  -s, --section       The section to run. Choices: -1, 0, 1, 2, 3.
  -c, --category      Define the category that the input sample belongs to. This will put the output files with the same category into the same container.
  -d, --isdata        Is data or not
```

### Check if the training inpnuts and reduced root-files are complete

The script `CheckArraysComplete.py` will check if all of the outputs from last step exist and create a text file which summarizes the missing files and can be used to recover the missing files.

```
python CheckArraysComplete.py
``` 

### Recover the missing files and bad files

After checking the completeness, if there are missing files, you can either run the recover script locally:

```
python RecoverMissingFiles.py -m
```

Or resubmit the jobs for the missing files:

```
python SubmitRecoverMissingFiles.py
```

### Start XGBoost analysis!

The whole ML task consists of training, applying the weights, plotting BDT distribution, optimizing the BDT boundaries for categorization, and calculating the number counting significances. The wrapper script `run_training.sh` will run everything. Please have a look!

#### Training a model

The training script `train_bdt.py` will train the model in four-fold, and transform the output scores such that the unweighted signal distribution is flat.

```
python train_bdt.py [-r REGION] [-f FOLD] [--VBF] [--ROC] [--save]

Usage:
  -r, --region        The region of interest. Choices: 'zero_jet', 'one_jet', or 'two_jet'.
  -f, --fold          Which fold to run. Default is -1 (run all folds)
  --VBF               To train VBF classifier
  --roc               To plot the ROC curve or not.
  --save              To save the model into HDF5 files, and the pickle files
```

#### Applying the weights

Applying the trained model (as well as the score transformation)  to the reduced ntuples to get BDT scores for each event can be done by doing:
```
python applyBDTWeight.py [-r REGION]
```
#### Optimizing the BDT boundaries

`categorization_1D.py` will take the Higgs classifier scores of the samples and optimized the boundaries that give the best combined significance. `submit_categorization_optimization_2D.py` will submit the jobs for optimizing the 2D scan of the BDT boundaries using both Higgs and VBF classifiers. Note that `submit_categorization_optimization_2D.py` can only be run after `categorization_1D.py` has been run at least once.

```
python categorization_1D.py [-r REGION] [-f NUMBER OF FOLDS] [-b NUMBER OF CATEGORIES] [-n NSCAN] [--floatB] [--minN minN]

Usage:
  -f, --fold          Number of folds of the categorization optimization. Default is 1.
  -b, --nbin          Number of BDT categories
  -n, --nscan         Number of scans. Default is 100
  --minN,             minN is the minimum number of events required in the mass window. The default is 5.
  --floatB            To float the last BDT boundary, which means to veto the lowest BDT score events
```

After the jobs for running 2D scan are done, one can use `getSignificance_2D.py` to extract the results.

```
python getSignificance_2D.py [-r REGION] [-b NUMBER OF ggF CATEGORIES] [-v NUMBER OF VBF CATEGORIES]
```
