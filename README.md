# Welcome to MonoHbbML

MonoHbbML is a package that can be used to study the XGBoost technique in Mono-H(bb) analysis.

## Setup

### First time setup on lxplus

First checkout the code:

```
git clone ssh://git@gitlab.cern.ch:7999/wisc_atlas/MonoHbbML.git [-b your_branch]
cd MonoHbbML
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
pip install theano keras h5py sklearn matplotlib tabulate xgboost
pip install --upgrade https://github.com/rootpy/root_numpy/zipball/master
```

If this is the first time you are using keras, you will want to change the backend to theano instead of the default tensorflow.
To do this, edit the appropriate line in `~/.keras/keras.json`.

### Normal setup on lxplus

After setting up the environment for the first time, you can return to this setup by doing `source setup_lxplus.sh`

## Scripts to run the tasks

### Check if the input files are complete

Before starting processing the arrays and all of the training analysis, it would be good to first check the completeness of the input ntuples. To do this, first setup rucio,

```
lsetup rucio
voms-proxy-init -voms atlas
```
and do `python checkinput.py`.

This script is no longer applicable as some files have been moved to other folders.

### Producing training inputs and reduced root-files which the trained model will be applied to

The core script is `process_arrays.py`. The script will apply the given preselection cuts to input samples and produce corresponding ntuples, numpy arrays files, and a text file containing the information of averaged weight of the sample. You can directly run it for any specified input files. Since there are too many samples in Mono-H(bb), it is more convinient and time-efficient to use `SubmitProcessArrays.py` which will find all of the input files and run the core script by submitting the condor jobs.


```
python SubmitProcessArrays.py
```
or
```
python process_arrays.py [-n INPUT_FILE_NAME] [-r REGION] [-v] [-t] [-c CATEGORY]

Note:
  -n, --name          The name of the input files containers. The mother folder ('inputs/') should not be included in the name. 
e.g. group.phys-exotics.mc16_13TeV.364222.Znunu_Sh221_PTV.0L0700a0.EXOT24.e5626_e5984_s3126_r9364_r9315_p3563_XAMPP/group.phys-exotics.14337810.XAMPP._000001.root
  -r, --region        The region of interest. 'm' for merged region. 'r_350_500' for the first resolved region. 'r_200_350' for the second resolved. 'r_150_200' for the last resolved.
  -v, --val           Do validation sample.
  -t, --train         Do training sample.
  -c, --category      Define the category that the input sample belongs to. This will put the output files of every sample that has the sample category into the same container.
```

### Check if the training inpnuts and reduced root-files are complete

The script `CheckArraysComplete.py` will check if all of the outputs from last step exist.

```
python CheckArraysComplete.py
``` 

### Recover the missing files and bad files

There are several stages where you will need to recover the files. One is when there are missing files. Others are when there are bad files that can't be read. Both cases can be addressed by running `RecoverMissingFiles.py`.

After checking the completeness, if there are missing files:

```
python RecoverMissingFiles.py -m
```

When running the training scripts later on, errors might occur due to some bad files that cannot be loaded. When this happens, do the following to recover the bad files:

```
python RecoverMissingFiles.py -b
```
Update: Now the later scripts are able to automatically call the recovering script to recover the bad files, so there might not be need anymore for doing the above command by hand.

### Start XGBoost analysis!

The whole ML task consists of training, applying the weights, plotting BDT distribution, optimizing the BDT boundaries for categorization, and calculating the number counting significances. `doTrain.py` will do everything for you, with assuming same signal model used for training and optimizing the categorization. Please note that `doTrain.py` only works for single trainged signal point and will only optimize the BDT boundaries using the same signal point. You have to run the following scripts separately to customze your training methods instead of simply running `doTrain.py`. 

```
python doTrain.py [-r REGION] [-m mZP mA]

```
#### Training a model

Now the codes in master branch work for multiple trained signals as well!
```
python train_bdt.py [-r REGION] [-m mZP_1 mA_1 mZP_2 mA_2 ...] [--save]
```
#### Applying the weights

Applying the trained model to the reduced ntuples to get BDT scores for each event can be done by doing:
```
python applyBDTWeight.py [-r REGION] [-m mZP_1 mA_1 mZP_2 mA_2 ...]
```
#### Plotting BDT distribution
`plot_score.py` will plot the BDT distribution for backgrounds and signal using their validation samples.
```
python plot_score.py [-r REGION] [-m Train_Model_mZP_1 Train_Model_mA_1 Train_Model_mA_2 ...] [-s Target_Signal_mZP_1 Target_Signal_mA_1 Target_signal_mZP_2 Target_signal_mA_2 ...]

Note:
  The trained signal points have to be in the same order as that you stated in the training.
```
#### Optimizing the BDT boundaries

Currently, one can only optimize the boundaries using one single signal point. The script in the master branch will get rid of leftmost bin (most background like) and output the best three BDT boundaries with which the combined significance of the target siganl is largest. Please checkout to branch full_cat if you want to also include the leftmost bin.
```
python categorization.py [-r REGION] [-m Train_Model_mZP_1 Train_Model_mA_1 Train_Model_mZP_2 Train_Model_mA_2 ...] [-s Target_Signal_mZP Target_Signal_mA]
```
#### Getting the significance using customized BDT boundaries

The number of the customized boundaries can be as many as you want. The boundaries have to be arranged in their ascending order.
```
python getSignificance.py [-r REGION] [-m Train_Model_mZP Train_Model_mA] [-b customized_boundary_1 customized_boundary_2 ...]
```

#### Plotting the significance and improvement
```
python plotSignificance.py [-r REGION] [-m Train_Model_mZP Train_Model_mA] [-s Target_Signal_mZP Target_Signal_mA] [-b customized_boundary_1 customized_boundary_2 ...]

Note:
  If you want to use the boundaries optimized on certain signal point, please use -s and specify the target signal point. 
  You can also plot the significances with your own customized BDT boundaries by using -b followed by the boundaries you prefer. 
  The number of the customized boundaries can be as many as you want. The boundaries have to be arranged in their ascending order.
``` 
#### Plotting distribution of each variable
```
python plotVariables.py [-r REGION] [-m mZP_1 mA_1 mZP_2 mA_2 ...]

Note:
  You can type as many masspoints of signals as you want.
```
#### Plotting feature-importance
```
python plot_features.py [-r REGION] [-m Train_Model_mZP_1 Train_Model_mA_1 Train_Model_mZP_2 Train_Model_mA_2 ...]
```
### Make a combined significance plot

You would need to modify to codes to use different models and BDT boundaries for each region.
```
python plotCombine.py
```
### Get the cut-based event yields
```
python getYields.py [-n name_to_specify_the_cut_methods]
```

### Get and plot baseline cut-based significances
```
python plotBaseline.py [-n name_of_the_cut_methods]
```
