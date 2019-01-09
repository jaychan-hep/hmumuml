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
import math
from argparse import ArgumentParser
import ROOT
import numpy as np
from array import array

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for Hmumu XGBoost analysis.")
    parser.add_argument('-n', '--name', action='store', required=True, help='File to process')
    parser.add_argument('-i', '--inputdir', action='store', default='inputs', help='Directory that contains ntuples.')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet','one_jet','zero_jet'],required=True, default='zero_jet', help='Region to process')
    parser.add_argument('-t', '--train', action='store_true', default=False, help='To process training samples')
    parser.add_argument('-v', '--val', action='store_true', default=False, help='To process validation samples')
    parser.add_argument('-c', '--category', action='store', required=True, help='The category of samples')
    parser.add_argument('-d', '--isdata', action='store_true', default=False, help='Real data') 
    return  parser.parse_args()

def process_arrays(args):

    DSID=args.name.split('.')[1]
    #container=args.name.split('/')[1]

    # Open input files
    f = ROOT.TFile.Open(args.inputdir+"/"+args.name)
    t = f.Get('DiMuonNtuple')

    # Define variables
    eventNumber = array('L', [0])

    nmuon = array('i',[0])
    Muons_Minv_MuMu = array('f', [0])
    Muons_PT_Lead = array('f', [0])
    Muons_PT_Sub = array('f', [0])
    Muons_Eta_Lead = array('f', [0])
    Muons_Eta_Sub = array('f', [0])
    Muons_Phi_Lead = array('f', [0])
    Muons_Phi_Sub = array('f', [0])

    Muons_MCPExpReso_Smear_Minv_MuMu_Sigma = array('f', [0])
    Muons_MCPExpReso_Smear_Lead = array('f', [0])
    Muons_MCPExpReso_Smear_Sub = array('f', [0])
    Muons_CosThetaStar = array('f', [0])

    Z_PT_FSR = array('f', [0])
    Z_Eta_FSR = array('f', [0])
    Z_Phi_FSR = array('f', [0])
    Z_Y_FSR = array('f', [0])

    njet = array('i',[0])
    njet_cen = array('i',[0])
    njet_fwd = array('i',[0])
    nbjet = array('i',[0])
    HT = array('f',[0])
    Jets_Pt = ROOT.vector('float')()
    Jets_Energy = ROOT.vector('float')()
    Jets_Eta = ROOT.vector('float')()
    Jets_Phi = ROOT.vector('float')()
    Jets_BTag = ROOT.vector('int')()
    Jets_DeltaEta_jj = array('f', [0])
    Jets_DeltaR_jj = array('f', [0])

    metFinalTrk = array('f', [0])
    metFinalTrkPhi = array('f', [0])
    Event_PT_MuMuj1 = array('f', [0])
    Event_PT_MuMuj2 = array('f', [0])
    Event_Y_MuMuj1 = array('f', [0])
    Event_Y_MuMuj2 = array('f', [0])
    Event_Centrality = array('f', [0])
    Event_MET_Sig = array('f', [0])

    if not args.isdata:
        Weight_Global = array('f', [0])
        Weight_Total = array('f', [0])


    #Set up the address of variables in branch
    t.SetBranchAddress('EventInfo_EventNumber', eventNumber)

    t.SetBranchAddress('Muons_Multip', nmuon)
    t.SetBranchAddress('Muons_Minv_MuMu', Muons_Minv_MuMu)
    t.SetBranchAddress('Muons_PT_Lead', Muons_PT_Lead)
    t.SetBranchAddress('Muons_PT_Sub', Muons_PT_Sub)
    t.SetBranchAddress('Muons_Eta_Lead', Muons_Eta_Lead)
    t.SetBranchAddress('Muons_Eta_Sub', Muons_Eta_Sub)
    t.SetBranchAddress('Muons_Phi_Lead', Muons_Phi_Lead)
    t.SetBranchAddress('Muons_Phi_Sub', Muons_Phi_Sub)

    t.SetBranchAddress('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma)
    t.SetBranchAddress('Muons_MCPExpReso_Smear_Lead', Muons_MCPExpReso_Smear_Lead)
    t.SetBranchAddress('Muons_MCPExpReso_Smear_Sub', Muons_MCPExpReso_Smear_Sub)
    t.SetBranchAddress('Muons_CosThetaStar', Muons_CosThetaStar)

    t.SetBranchAddress('Z_PT_FSR', Z_PT_FSR)
    t.SetBranchAddress('Z_Eta_FSR', Z_Eta_FSR)
    t.SetBranchAddress('Z_Phi_FSR', Z_Phi_FSR)
    t.SetBranchAddress('Z_Y_FSR', Z_Y_FSR)

    #t.SetBranchAddress('Jets_Multip', njet)
    t.SetBranchAddress('Jets_PT', Jets_Pt)
    t.SetBranchAddress('Jets_E', Jets_Energy)
    t.SetBranchAddress('Jets_Eta', Jets_Eta)
    t.SetBranchAddress('Jets_Phi', Jets_Phi)
    t.SetBranchAddress('Jets_LowestPassedBTagOP', Jets_BTag)
    t.SetBranchAddress('Jets_DeltaEta_jj', Jets_DeltaEta_jj)
    t.SetBranchAddress('Jets_DeltaR_jj', Jets_DeltaR_jj)

    t.SetBranchAddress('Event_MET', metFinalTrk)
    t.SetBranchAddress('Event_MET_Phi', metFinalTrkPhi)
    t.SetBranchAddress('Event_PT_MuMuj1', Event_PT_MuMuj1)
    t.SetBranchAddress('Event_PT_MuMuj2', Event_PT_MuMuj2)
    t.SetBranchAddress('Event_Y_MuMuj1', Event_Y_MuMuj1)
    t.SetBranchAddress('Event_Y_MuMuj2', Event_Y_MuMuj1)
    t.SetBranchAddress('Event_Centrality', Event_Centrality)
    t.SetBranchAddress('Event_MET_Sig', Event_MET_Sig)

    if not args.isdata:
        t.SetBranchAddress('GlobalWeight', Weight_Global)
        t.SetBranchAddress('TotalWeight', Weight_Total)

    # to run faster
    t.SetBranchStatus("*",0)

    t.SetBranchStatus('EventInfo_EventNumber', 1)

    t.SetBranchStatus('Muons_Multip', 1)
    t.SetBranchStatus('Muons_Minv_MuMu', 1)
    t.SetBranchStatus('Muons_PT_Lead', 1)
    t.SetBranchStatus('Muons_PT_Sub', 1)
    t.SetBranchStatus('Muons_Eta_Lead', 1)
    t.SetBranchStatus('Muons_Eta_Sub', 1)
    t.SetBranchStatus('Muons_Phi_Lead', 1)
    t.SetBranchStatus('Muons_Phi_Sub', 1)

    t.SetBranchStatus('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', 1)
    t.SetBranchStatus('Muons_MCPExpReso_Smear_Lead', 1)
    t.SetBranchStatus('Muons_MCPExpReso_Smear_Sub', 1)
    t.SetBranchStatus('Muons_CosThetaStar', 0)
    t.SetBranchStatus('Z_PT_FSR', 1)
    t.SetBranchStatus('Z_Eta_FSR', 1)
    t.SetBranchStatus('Z_Phi_FSR', 1)
    t.SetBranchStatus('Z_Y_FSR', 1)

    #t.SetBranchStatus('Jets_Multip', 1)
    t.SetBranchStatus('Jets_PT', 1)
    t.SetBranchStatus('Jets_E', 1)
    t.SetBranchStatus('Jets_Eta', 1)
    t.SetBranchStatus('Jets_Phi', 1)
    t.SetBranchStatus('Jets_LowestPassedBTagOP', 1)
    t.SetBranchStatus('Jets_DeltaEta_jj', 1)

    t.SetBranchStatus('Event_MET', 1)
    t.SetBranchStatus('Event_MET_Phi', 1)
    t.SetBranchStatus('Event_PT_MuMuj1', 1)
    t.SetBranchStatus('Event_PT_MuMuj2', 0)
    t.SetBranchStatus('Event_Y_MuMuj1', 1)
    t.SetBranchStatus('Event_Y_MuMuj2', 1)
    t.SetBranchStatus('Event_Centrality', 0)
    t.SetBranchStatus('Event_MET_Sig', 0)

    if not args.isdata:
        t.SetBranchStatus('GlobalWeight', 1)
        t.SetBranchStatus('TotalWeight', 1)

    if not args.train and not args.val:
        # Create new tree
        outtree = ROOT.TTree('AppInput', 'AppInput')
        outtree.SetDirectory(0)

        # Create branches for tree
        weight=array('d',[0])
        m_mumu=array('f',[0])
        pt_mumu=array('f',[0])

        outtree.Branch('eventNumber', eventNumber, 'eventNumber/l')
        outtree.Branch('njet', njet, 'njet/I')
        outtree.Branch('nbjet', nbjet, 'nbjet/I')
        outtree.Branch('metFinalTrk', metFinalTrk, 'metFinalTrk/F')
        outtree.Branch('m_mumu', m_mumu, 'm_mumu/F')
        outtree.Branch('pt_mumu', pt_mumu, 'pt_mumu/F')
        # if not args.isdata:
        outtree.Branch('weight', weight, 'weight/D')

    met = []
    jets = []
    muons = []
    dijet = []
    dimuon = []
    extras = []
    if args.train or args.val: # and not args.isdata:
        wt = []

    print "Running %d events" % (t.GetEntries())
    
    event_yield=0
    #cutflow_check=0
    for i in range(t.GetEntries()):
        if (i % 10000 == 0): print('Event %d/%d' % (i, t.GetEntries()))
        t.GetEntry(i)
        
        #if int(DSID) <= 364197 and int(DSID) >= 364100:
            #if not abs(mcEventWeight[0]) < 100: continue

        # place holder 
        if 1:
            # preselection
            if not (nmuon[0] == 2): continue
            if not (Muons_Minv_MuMu[0] >= 105 and Muons_Minv_MuMu[0] <= 180): continue
            if (args.train or args.val):
                if args.isdata:
                    if (Muons_Minv_MuMu[0] >= 120 and Muons_Minv_MuMu[0] <= 130): continue
                else:
                    if not (Muons_Minv_MuMu[0] >= 120 and Muons_Minv_MuMu[0] <= 130): continue

            # separate events into training and validation samples
            if (args.train and eventNumber[0]%2 != 0): continue
            if (args.val and eventNumber[0]%2 != 1): continue

            # mets
            eM = []
            eM.append([metFinalTrk[0], metFinalTrkPhi[0]])

            # Jets
            ej = []
            njet[0] = 0
            njet_cen[0] = 0
            njet_fwd[0] = 0
            nbjet[0] = 0
            HT[0] = 0
            for j in range(0, len(Jets_Pt)):
                njet[0] = njet[0] + 1
                HT[0] = HT[0] + Jets_Pt[j]
                if (abs(Jets_Eta[j]) < 2.5): njet_cen[0] = njet_cen[0] + 1
                else: njet_fwd[0] = njet_fwd[0] + 1
                if (Jets_BTag[j] == 60): 
                    nbjet[0] = nbjet[0] + 1
                    ej.append([Jets_Pt[j], Jets_Eta[j], Jets_Phi[j], 0, 1])
                else: 
                    ej.append([Jets_Pt[j], Jets_Eta[j], Jets_Phi[j], 0, 0])
            ej = sorted(ej, key=lambda x: x[0], reverse=True)
            if (len(ej) > 2):
                for j in range(len(ej) - 2): ej.pop()
            elif (len(ej) < 2):
                for j in range(2 - len(ej)): ej.append([-999] * 5)

            Jets_1 = ROOT.TLorentzVector()
            if(njet[0] >= 1): Jets_1.SetPtEtaPhiM(ej[0][0], ej[0][1], ej[0][2], 0)

            Jets_2 = ROOT.TLorentzVector()
            if(njet[0] >= 2): Jets_2.SetPtEtaPhiM(ej[1][0], ej[1][1], ej[1][2], 0)

            j1j2 = ROOT.TLorentzVector()
            if(njet[0] >= 2): j1j2 = Jets_1+Jets_2

            # muons
            Muons_Lead = ROOT.TLorentzVector()
            Muons_Lead.SetPtEtaPhiM(Muons_PT_Lead[0], Muons_Eta_Lead[0], Muons_Phi_Lead[0], 0.106)

            Muons_Sub = ROOT.TLorentzVector()
            Muons_Sub.SetPtEtaPhiM(Muons_PT_Sub[0], Muons_Eta_Sub[0], Muons_Phi_Sub[0], 0.106)

            mumu = ROOT.TLorentzVector()
            mumu = Muons_Lead+Muons_Sub

            em = []
            em.append([Muons_PT_Lead[0]/Muons_Minv_MuMu[0], Muons_Eta_Lead[0], Muons_Phi_Lead[0], 0])
            em.append([Muons_PT_Sub[0]/Muons_Minv_MuMu[0], Muons_Eta_Sub[0], Muons_Phi_Sub[0], 0])

            if not args.train and not args.val:
                m_mumu[0] = Muons_Minv_MuMu[0]
                pt_mumu[0] = mumu.Pt()

            # dijet
            dj = []
            if (njet[0] >= 2): dj.append([j1j2.Pt(), j1j2.Eta(), j1j2.Phi(), j1j2.M(), 0, 0])
            else: dj.append([-999] * 6)

            # dimuon
            dm = []
            dm.append([Z_PT_FSR[0]/Muons_Minv_MuMu[0], Z_Eta_FSR[0], Z_Phi_FSR[0], Z_Y_FSR[0], 0])

            # extras
            ee = []
            ee.append([njet_cen[0], njet_fwd[0], nbjet[0], HT[0], Event_PT_MuMuj1[0], Event_Y_MuMuj1[0], 0])

            # cuts
            if not (nbjet[0] == 0 and metFinalTrk[0] < 80): continue
            if args.region == 'two_jet':
                if not (njet[0] >= 2): continue
            if args.region == 'one_jet':
                if not (njet[0] == 1): continue
            if args.region == 'zero_jet':
                if not (njet[0] == 0): continue
            event_yield+=1

            met.append(eM)
            jets.append(ej)
            muons.append(em)
            dijet.append(dj)
            dimuon.append(dm)
            extras.append(ee)
            
            # weight
            if not args.isdata:
                if not args.train and not args.val:
                    weight[0] = Weight_Global[0]
                else:
                    wt.append(Weight_Global[0])
            else:
                if not args.train and not args.val:
                    weight[0] = 1
                else:
                    wt.append(1)

        if not args.train and not args.val:
            outtree.Fill()
    

    print "INFO: The event yield after preselection is: ", event_yield

    # If the arrays are not empty, then convert them into numpy arrays and save them

    #output_container = args.name.split('/')[0]
    category=args.category     
    if not os.path.isdir('arrays'):
        print 'INFO: Creating output folder: \"arrays\"'
        os.makedirs("arrays")
    if not os.path.isdir('arrays/'+category):
        print 'INFO: Creating new arrays container: '+ category
        os.makedirs("arrays/"+category)
        
  
    array_basename=DSID
    if 'MC16A' in args.name or 'mc16a' in args.name:
        array_basename='mc16a_'+array_basename
    elif 'MC16D' in args.name or 'mc16d' in args.name:
        array_basename='mc16d_'+array_basename
    elif 'data15' in args.name:
        array_basename='data15_'+array_basename
    elif 'data16' in args.name:
        array_basename='data16_'+array_basename
    elif 'data17' in args.name:
        array_basename='data17_'+array_basename
    else:
        print 'WARNING: File does not match with a tipical file name!!'

    if args.region == 'two_jet':
        region='two_jet'
    elif args.region == 'one_jet':
        region='one_jet'
    else:
        region='zero_jet'

    if event_yield==0:
        print "WARNING: No event left after the preselections. Will not save the results."
        with open('arrays/%s/%s%s%s_%s.txt' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region),"w") as txt:
            txt.write('%s%s%s_%s.npy\n' % (array_basename, '_train' if args.train else '', '_val' if args.val else '', region))
    else:
        if args.train or args.val: # and not args.isdata:
            wt = np.array(wt)
            print 'Saving numpy arrays of weight to "arrays/%s"...' % (category)
            np.save('arrays/%s/%s%s%s_%s_weight.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), wt)

        met = np.array(met)
        print 'Saving numpy arrays of MET to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s%s%s_%s_met.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), met)

        jets = np.array(jets)
        print 'Saving numpy arrays of jets to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s%s%s_%s_jets.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), jets)

        muons = np.array(muons)
        print 'Saving numpy arrays of muons to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s%s%s_%s_muons.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), muons)

        dijet = np.array(dijet)
        print 'Saving numpy arrays of dijet to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s%s%s_%s_dijet.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), dijet)

        dimuon = np.array(dimuon)
        print 'Saving numpy arrays of dimuon to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s%s%s_%s_dimuon.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), dimuon)

        extras = np.array(extras)
        print 'Saving numpy arrays of extras to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s%s%s_%s_extras.npy' % (category, array_basename, '_train' if args.train else '', '_val' if args.val else '', region), extras)

        if not args.train and not args.val:
            if not os.path.isdir('AppInputs'):
                print 'INFO: Creating new folder: \"AppInputs\"'
                os.makedirs("AppInputs")
            if not os.path.isdir('AppInputs/'+category):
                print 'INFO: Creating new AppInputs container: '+ category
                os.makedirs("AppInputs/"+category)

            print 'Saving AppInputs to \"AppInputs/%s\"...' % (category)
            outfile=ROOT.TFile('AppInputs/%s/%s_%s.root' %(category, array_basename, region), 'RECREATE')
            outtree.Write()
            outfile.Close()

    return


def main():
    """Run numpy arrays production for Hmumu XGBoost analysis."""
    args=getArgs()
    if args.train and args.val: 
        print 'ERROR: Please select only one running section!! (either --val or --train)'
        print 'Exiting...'
        quit()
    if not (args.region == 'two_jet' or args.region == 'one_jet' or args.region == 'zero_jet'):
        print 'ERROR: The region is not defined correctly. Please select from \'two_jet\', \'one_jet\', \'zero_jet\'.'
        print 'Exiting...'
        quit()
    if not os.path.isfile("inputs/"+args.name):
        print "ERROR: The input file %s doesn\'t exist!!" % ("inputs/"+args.name)
        print 'Exiting...'
        quit()
    if args.train:
        section='training'
    elif args.val:
        section='validation'
    else:
        section='general'
    print "========================================================="
    print "Processing %s" % args.name
    print "========================================================="
    print ""
    print "INFO: Making %s samples" % section
    if args.region =='two_jet': print "INFO: Preselection for the two jet region."
    elif args.region =='one_jet': print "INFO: Preselection for the one jet region."
    else: print "INFO: Preselection for the zero jet region."

    process_arrays(args)
    
    print ""
    print "========================================================="
    print "Finishing the process"
    print "=========================================================" 
    return

if __name__ == '__main__':
    main()
