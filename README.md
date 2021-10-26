# Welcome to HmumuML

HmumuML is a package that can be used to study the XGBoost BDT categorization in Hmumu analysis.

## Setup

This project requires to use python3 and either conda or virtual environment. If you want to use a conda environment, install Anaconda or Miniconda before setting up.

### First time setup on lxplus for a conda environment

First checkout the code:

```
git clone ssh://git@gitlab.cern.ch:7999/HZZ/HZZSoftware/Hmumu/hmumuml.git [-b your_branch]
cd HmumuML
```

Then,

```
source scripts/install.sh
```
### First time setup on lxplus for a virtual environment
After cloning the repository, under HmumuML:
```
virtualenv hmumuml --python=python3
source hmumuml/bin/activate
pip install -r requirements.txt
```

### Normal setup on lxplus

After setting up the environment for the first time, you can return to this setup by doing `source scripts/setup.sh`

## Scripts to run the tasks

### Prepare the training inputs (skimmed ntuples)

The core script is `skim_ntuples.py`. The script will apply the given skimming cuts to input samples and produce corresponding skimmed ntuples. You can run it locally for any specified input files. In H->mumu, it is more convenient and time-efficient to use `submit_skim_ntuples.py` which will find all of the input files specified in `data/inputs_config.json` and run the core script by submitting the condor jobs.

- The script is hard-coded currently, meaning one needs to directly modify the core script to change the variable calculation and the skimming cuts.
- The output files (skimmed ntuples) will be saved to the folder named `skimmed_ntuples` by default.
- `submit_skim_ntuples.py` will only submit the jobs for those files that don't exist in the output folder.

```
python scripts/submit_skim_ntuples.py
```
or
```
python scripts/skim_ntuples.py [-i input_file_path] [-o output_file_path]
```

### Check the completeness of the jobs or recover the incomplete jobs

The script `run_skim_ntuples.py` will check if all of the outputs from last step exist and if the resulting number of events is correct.

```
python scripts/run_skim_ntuples.py [-c]

usage:
  -c   check the completeness only (without recovering the incomplete jobs)
  -s   skip checking the resulting number of events
``` 

Please be sure to do `python run_skim_ntuples.py -c` at least once before starting the BDT training exercise.


### Start XGBoost analysis!

The whole ML task consists of training, applying the weights, optimizing the BDT boundaries for categorization, and calculating the number counting significances. The wrapper script `run_all.sh` will run everything. Please have a look!

#### Training a model

The training script `train_bdt.py` will train the model in four-fold, and transform the output scores such that the unweighted signal distribution is flat. The detailed settings, including the preselections, training variables, hyperparameters, etc, are specified in the config file `data/training_config.json`.

```
python scripts/train_bdt.py [-r TRAINED_MODEL] [-f FOLD] [--save]

Usage:
  -r, --region        The model to be trained. Choices: 'zero_jet', 'one_jet', 'two_jet' or 'VBF'.
  -f, --fold          Which fold to run. Default is -1 (run all folds)
  --save              To save the model into HDF5 files, and the pickle files
```

#### Applying the weights

Applying the trained model (as well as the score transformation) to the skimmed ntuples to get BDT scores for each event can be done by doing:
```
python scripts/apply_bdt.py [-r REGION]
```
The script will take the settings specified in the training config file `data/training_config.json` and the applying config file `data/apply_config.json`.

#### Optimizing the BDT boundaries

`categorization_1D.py` will take the Higgs classifier scores of the samples and optimize the boundaries that give the best combined significance. `categorization_2D.py`, on the other hand, takes both the Higgs classifier scores and the VBF classifier scores of the samples and optimizes the 2D boundaries that give the best combined significance.

```
python categorization_1D.py [-r REGION] [-f NUMBER OF FOLDS] [-b NUMBER OF CATEGORIES] [-n NSCAN] [--floatB] [--minN minN] [--skip]

Usage:
  -f, --fold          Number of folds of the categorization optimization. Default is 1.
  -b, --nbin          Number of BDT categories
  -n, --nscan         Number of scans. Default is 100
  --minN,             minN is the minimum number of events required in the mass window. The default is 5.
  --floatB            To float the last BDT boundary, which means to veto the lowest BDT score events
  --skip              To skip the hadd step (if you have already merged signal and background samples)
```

```
python categorization_2D.py [-r REGION] [-f NUMBER OF FOLDS] [-b NUMBER OF CATEGORIES] [-b NUMBER OF ggF CATEGORIES] [-n NSCAN] [--floatB] [--minN minN] [--skip]

Usage:
  -f, --fold          Number of folds of the categorization optimization. Default is 1.
  -b, --nbin          Number of BDT categories
  -n, --nscan         Number of scans. Default is 100
  --minN,             minN is the minimum number of events required in the mass window. The default is 5.
  --floatB            To float the last BDT boundary, which means to veto the lowest BDT score events
  --skip              To skip the hadd step (if you have already merged signal and background samples)
```
