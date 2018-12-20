#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.12.2018
#
#  Modified from ttHyy XGBosst anaylsis codes
#
#
#
#
import os
from argparse import ArgumentParser
import math
import ROOT
from root_numpy import root2array
from array import array
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
from keras import models
import xgboost as xgb


def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet'], default='zero_jet', help='Region to process')
    return parser.parse_args()


def apply_weight(train_sig,catagories,region,isdata=False):

    model = xgb.Booster()
    model.load_model('models/%s%s.h5'%(train_sig, region))

    model_r = xgb.Booster()
    model_r.load_model('models/%s%s_r.h5'%(train_sig, region))

    # loading
    badfile=False
    for category in catagories:
        print 'Reading from %s samples...' %(category)
        AppInputsdir="AppInputs/"+category
        npydir='arrays/'+category
        for appinput in os.listdir(AppInputsdir):
            if not appinput.endswith('_%s.root' % (region)):
                continue
            # if '_train_' in appinput or '_val_' in appinput: continue
            try:
                jets=np.load(npydir+'/'+appinput.replace('.root','_jets.npy'))
                if [jets.shape[1],jets.shape[2]]!=[2,5]: 
                    print 'ERROR: Dimension of jets npy is not correct!!'
                    error
                met=np.load(npydir+'/'+appinput.replace('.root','_met.npy'))
                if [met.shape[1],met.shape[2]]!=[1,2]:
                    print 'ERROR: Dimension of met npy is not correct!!'
                    error
                muons=np.load(npydir+'/'+appinput.replace('.root','_muons.npy'))
                if [muons.shape[1],muons.shape[2]]!=[2,4]:
                    print 'ERROR: Dimension of muons npy is not correct!!'
                    error
                dijet=np.load(npydir+'/'+appinput.replace('.root','_dijet.npy'))
                if [dijet.shape[1],dijet.shape[2]]!=[1,6]:
                    print 'ERROR: Dimension of dijet npy is not correct!!'
                    error
                dimuon=np.load(npydir+'/'+appinput.replace('.root','_dimuon.npy'))
                if [dimuon.shape[1],dimuon.shape[2]]!=[1,5]:
                    print 'ERROR: Dimension of dimuon npy is not correct!!'
                    error
                extras=np.load(npydir+'/'+appinput.replace('.root','_extras.npy'))
                if [extras.shape[1],extras.shape[2]]!=[1,6]:
                    print 'ERROR: Dimension of extras npy is not correct!!'
                    error
                infile = ROOT.TFile.Open(AppInputsdir+'/'+appinput)
                intree = infile.Get('AppInput')
        
                eventNumber = array('L', [0])
                m_mumu = array('f', [0])
                njet = array('I', [0])
                nbjet = array('I', [0])
                metFinalTrk = array('f', [0])
                if not isdata:
                    weight=array('d',[0])
                intree.SetBranchAddress('eventNumber', eventNumber)
                intree.SetBranchAddress('m_mumu', m_mumu)
                intree.SetBranchAddress('njet', njet)
                intree.SetBranchAddress('nbjet', nbjet)
                intree.SetBranchAddress('metFinalTrk', metFinalTrk)
                if not isdata:
                    intree.SetBranchAddress('weight', weight)

                # output tree
                bdt_score = array('f', [0])
                outtree = ROOT.TTree('test', 'test')
                outtree.SetDirectory(0)
                outtree.Branch('eventNumber', eventNumber,'eventNumber/l')
                outtree.Branch('m_mumu', m_mumu,'m_mumu/F')
                outtree.Branch('njet', njet, 'njet/i')
                outtree.Branch('nbjet', nbjet, 'nbjet/i')
                outtree.Branch('metFinalTrk', metFinalTrk,'metFinalTrk/F')
                if not isdata:
                    outtree.Branch('weight', weight, 'weight/D')
                outtree.Branch('bdt_score', bdt_score, 'bdt_score/F')
    
                # get BDT scores for given arrays
                #print 'Getting BDT scores from the trained model...'
                events=[]
                events.append(jets.reshape((jets.shape[0], -1)))
                events.append(met.reshape((met.shape[0], -1)))
                events.append(muons.reshape((muons.shape[0], -1)))
                events.append(dijet.reshape((dijet.shape[0], -1)))
                events.append(dimuon.reshape((dimuon.shape[0], -1)))
                events.append(extras.reshape((extras.shape[0], -1)))
                events = np.concatenate(events, axis=1)
                dEvents = xgb.DMatrix(events)
 
                score = model.predict(dEvents)
                score_r = model_r.predict(dEvents)

                if len(score) != intree.GetEntries():
                    print "ERROR: The number of events doesn't match!!!"
                    quit()

                #print 'writing bdt scores into new rootfiles...'
                for i in range(intree.GetEntries()):
                    intree.GetEntry(i)
                    if (eventNumber[0]%2==1): 
                        bdt_score[0] = 1.0 - score[i]
                    else:
                        bdt_score[0] = 1.0 - score_r[i]
                    outtree.Fill()

                if not os.path.isdir('outputs'):
                    print 'INFO: Creating new folder: \"outputs\"'
                    os.makedirs("outputs")
                if not os.path.isdir('outputs/model_%s%s' %(train_sig,region)):
                    print 'INFO: Creating model container:  model_%s%s' %(train_sig,region)
                    os.makedirs('outputs/model_%s%s' %(train_sig,region))
                if not os.path.isdir('outputs/model_%s%s/%s' %(train_sig,region,category)):
                    print 'INFO: Creating in the model container a new category:  %s' %(category)
                    os.makedirs('outputs/model_%s%s/%s' %(train_sig,region,category))


                outfile = ROOT.TFile('outputs/model_%s%s/%s/%s' % (train_sig,region,category,appinput), 'RECREATE')
                outtree.Write()
                outfile.Close()

                infile.Close()
            except:
                badfile=True
                print('ERROR: Unable to read \"'+ AppInputsdir+'/'+appinput+'\"')
                #stop=0
                #for container in os.listdir('inputs'):
                #    if stop: break
                #    if category in container:
                #        for ntuple in os.listdir('inputs/'+container):
                #            if appinput.split('_')[1] in ntuple and appinput.split('_')[2] in ntuple:
                #                print 'Original file: \"%s/%s\"'%(container,ntuple)
                #                try:
                #                    if txt: pass
                #                except:
                #                    txt = open('arrays/badsamplelist.txt','w')
                #                txt.write('%s %s %s %s %s badsample\n'%(container, ntuple, 'none', region, category))
                #                stop=1
                #                break

        if badfile:
            return False


        print 'Combining all the root-files in the category:  ', category
        if os.path.isdir('outputs/model_%s%s/%s'%(train_sig,region,category)):
            mc16a=False
            mc16d=False
            for out in os.listdir('outputs/model_%s%s/%s'%(train_sig,region,category)):
                if 'mc16a' in out: mc16a=True
                if 'mc16d' in out: mc16d=True
                if mc16a and mc16d: break
            #if mc16a:
            #    os.system("hadd -f outputs/model_%s%s/%s_mc16a.root outputs/model_%s%s/%s/mc16a*.root" %(train_sig,region,category,train_sig,region,category))
            #if mc16d:
            #    os.system("hadd -f outputs/model_%s%s/%s_mc16d.root outputs/model_%s%s/%s/mc16d*.root" %(train_sig,region,category,train_sig,region,category))
            os.system("hadd -f outputs/model_%s%s/%s.root outputs/model_%s%s/%s/*.root" %(train_sig,region,category,train_sig,region,category))

    return True



def main():
    args=getArgs()
    

    sigs = ['ggF','VBF']#,'VH','ttH']
    train_sig=''
    bkgs = ['data']

    if args.region == 'two_jet':
        region = 'two_jet'
    elif args.region == 'one_jet':
        region = 'one_jet'
    else:
        region = 'zero_jet'

    print '=============================================================================='
    print 'INFO:  Using the model: %s%s.h'%(train_sig,region) 
    print '------------------------------------------------------------------------------'

    while not apply_weight(train_sig,sigs,region):
        os.system('python RecoverMissingFiles.py -b')
    while not apply_weight(train_sig,bkgs,region):
        os.system('python RecoverMissingFiles.py -b')

    print '------------------------------------------------------------------------------'
    print 'Finishing the process'

    return

if __name__ == '__main__':
    main()
