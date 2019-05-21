#!/bin/bash
#########################################################################
#                                                                       #
#    The example wrapper for the training and categorization task.      #
#                                                                       #
#    Please contact ~jay.chan@cern.ch~ in case there is any issue.      #
#                                                                       #
#########################################################################

############################
#  Training the BDT models
############################
python train_bdt.py -r zero_jet --save
python train_bdt.py -r one_jet --save
python train_bdt.py -r two_jet --save
python train_bdt.py -r two_jet --VBF --save

###########################################
#  Applying the BDT models to all samples
###########################################
python applyBDTWeight.py -r zero_jet
python applyBDTWeight.py -r one_jet
python applyBDTWeight.py -r two_jet

###########################################################
#  Optimizing the BDT boundaries for zero-jet and two-jet
###########################################################
python categorization_1D.py -r zero_jet -b 3
python categorization_1D.py -r one_jet -b 3

##############################################
#  Optimizing the BDT boundaries for two-jet
##############################################
python categorization_1D.py -r two_jet -b 3
python submit_categorization_optimization_2D.py

#################################################################################################################################
# After all the jobs are done, one needs to do `python getSignificance_2D.py -r two_jet -v 3 -b 3` to get the result of two-jet
################################################################################################################################
