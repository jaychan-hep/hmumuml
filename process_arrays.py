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
import time
import sys

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for Hmumu XGBoost analysis.")
    parser.add_argument('-n', '--name', action='store', required=True, help='File to process')
    parser.add_argument('-i', '--inputdir', action='store', default='inputs', help='Directory that contains ntuples.')
    parser.add_argument('-r', '--region', action='store', choices=['all_jet','two_jet','one_jet','zero_jet'],default='all_jet', help='Region to process')
    parser.add_argument('-s', '--section', type=int, choices=[-1, 0, 1, 2, 3], required=True, help='The section to process')
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
    Muons_Minv_MuMu_Fsr = array('f', [0])

    Muons_MCPExpReso_Smear_Minv_MuMu_Sigma = array('f', [0])
    Muons_Minv_MuMu_Sigma = array('f', [0])
    Muons_CosThetaStar = array('f', [0])

    Z_PT_FSR = array('f', [0])
    Z_Y_FSR = array('f', [0])
    Z_Phi_FSR = array('f', [0])

    Jets_PT_Lead = array('f', [0])
    Jets_PT_Sub = array('f', [0])
    Jets_Eta_Lead = array('f', [0])
    Jets_Eta_Sub = array('f', [0])
    Jets_Phi_Lead = array('f', [0])
    Jets_Phi_Sub = array('f', [0])
    DeltaPhi_mumuj1 = array('f', [0])
    DeltaPhi_mumuj2 = array('f', [0])
    Jets_PT_jj = array('f', [0])
    Jets_Eta_jj = array('f', [0])
    Jets_Y_jj = array('f', [0])
    Jets_Phi_jj = array('f', [0])
    DeltaPhi_mumujj = array('f', [0])
    Jets_Minv_jj = array('f', [0])

    njet = array('i',[0])

    metFinalTrk = array('f', [0])

    Event_HasBJet = array('i', [0])

    Weight_Global = array('f', [0])

    #Set up the address of variables in branch
    t.SetBranchAddress('EventInfo_EventNumber', eventNumber)

    t.SetBranchAddress('Muons_Multip', nmuon)
    t.SetBranchAddress('Muons_Minv_MuMu_Fsr', Muons_Minv_MuMu_Fsr)

    t.SetBranchAddress('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma)
    t.SetBranchAddress('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma)
    t.SetBranchAddress('Muons_CosThetaStar', Muons_CosThetaStar)

    t.SetBranchAddress('Z_PT_FSR', Z_PT_FSR)
    t.SetBranchAddress('Z_Y_FSR', Z_Y_FSR)
    t.SetBranchAddress('Z_Phi_FSR', Z_Phi_FSR)

    t.SetBranchAddress('Jets_PT_Lead', Jets_PT_Lead)
    t.SetBranchAddress('Jets_PT_Sub', Jets_PT_Sub)
    t.SetBranchAddress('Jets_Eta_Lead', Jets_Eta_Lead)
    t.SetBranchAddress('Jets_Eta_Sub', Jets_Eta_Sub)
    t.SetBranchAddress('Jets_Phi_Lead', Jets_Phi_Lead)
    t.SetBranchAddress('Jets_Phi_Sub', Jets_Phi_Sub)
    t.SetBranchAddress('Jets_PT_jj', Jets_PT_jj)
    t.SetBranchAddress('Jets_Eta_jj', Jets_Eta_jj)
    t.SetBranchAddress('Jets_Phi_jj', Jets_Phi_jj)
    t.SetBranchAddress('Jets_Minv_jj', Jets_Minv_jj)

    t.SetBranchAddress('Jets_jetMultip', njet)

    t.SetBranchAddress('Event_HasBJet', Event_HasBJet)
    t.SetBranchAddress('Event_MET', metFinalTrk)

    t.SetBranchAddress('GlobalWeight', Weight_Global)

    if args.section == -1:
        # Create new tree
        outtree = ROOT.TTree('AppInput', 'AppInput')
        outtree.SetDirectory(0)

        # Create branches for tree
        weight=array('d',[0])
        m_mumu=array('f',[0])
        pt_mumu=array('f',[0])

        outtree.Branch('eventNumber', eventNumber, 'eventNumber/l')
        outtree.Branch('njet', njet, 'njet/I')

        outtree.Branch('metFinalTrk', metFinalTrk, 'metFinalTrk/F')

        outtree.Branch('Z_PT_FSR', Z_PT_FSR, 'Z_PT_FSR/F')
        outtree.Branch('Z_Y_FSR', Z_Y_FSR, 'Z_Y_FSR/F')
        outtree.Branch('Z_Phi_FSR', Z_Phi_FSR, 'Z_Phi_FSR/F')

        outtree.Branch('Jets_PT_Lead', Jets_PT_Lead, 'Jets_PT_Lead/F')
        outtree.Branch('Jets_PT_Sub', Jets_PT_Sub, 'Jets_PT_Sub/F')
        outtree.Branch('Jets_Eta_Lead', Jets_Eta_Lead, 'Jets_Eta_Lead/F')
        outtree.Branch('Jets_Eta_Sub', Jets_Eta_Sub, 'Jets_Eta_Sub/F')
        outtree.Branch('DeltaPhi_mumuj1', DeltaPhi_mumuj1, 'DeltaPhi_mumuj1/F')
        outtree.Branch('DeltaPhi_mumuj2', DeltaPhi_mumuj2, 'DeltaPhi_mumuj2/F')

        outtree.Branch('Jets_PT_jj', Jets_PT_jj, 'Jets_PT_jj/F')
        outtree.Branch('Jets_Y_jj', Jets_Y_jj, 'Jets_Y_jj/F')
        outtree.Branch('DeltaPhi_mumujj', DeltaPhi_mumujj, 'DeltaPhi_mumujj/F')
        outtree.Branch('Jets_Minv_jj', Jets_Minv_jj, 'Jets_Minv_jj/F')

        outtree.Branch('Muons_CosThetaStar', Muons_CosThetaStar, 'Muons_CosThetaStar/F')

        outtree.Branch('Muons_MCPExpReso_Smear_Minv_MuMu_Sigma', Muons_MCPExpReso_Smear_Minv_MuMu_Sigma, 'Muons_MCPExpReso_Smear_Minv_MuMu_Sigma/F')
        outtree.Branch('Muons_Minv_MuMu_Sigma', Muons_Minv_MuMu_Sigma, 'Muons_Minv_MuMu_Sigma/F')

        outtree.Branch('m_mumu', Muons_Minv_MuMu_Fsr, 'm_mumu/F')
        outtree.Branch('weight', weight, 'weight/D')

    met = []
    jets = []
    dijet = []
    dimuon = []
    extras = []
    if args.section != -1: # and not args.isdata:
        wt = []

    print "INFO: Running %d events" % (t.GetEntries())
    
    event_yield=0
    #cutflow_check=0

    start_time = time.time()

    for i in range(t.GetEntries()):

        if (i % 1000 == 0): 
            elapsed_time = time.time() - start_time
            tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
            print '%s |%s%s| %d/100'%(tt, '>'*(i*100/t.GetEntries()), ' '*(100-i*100/t.GetEntries()), i*100/t.GetEntries())
            sys.stdout.write("\033[F") # Cursor up one line

        t.GetEntry(i)
        
        # place holder 
        if 1:
            # preselection
            if args.region == 'two_jet':
                if not (njet[0] >= 2): continue
            if args.region == 'one_jet':
                if not (njet[0] == 1): continue
            if args.region == 'zero_jet':
                if not (njet[0] == 0): continue

            if not (nmuon[0] == 2): continue
            if not (Event_HasBJet[0] == 0): continue
            #if not (args.region == 'two_jet' or metFinalTrk[0] <= 80): continue
            if not (Muons_Minv_MuMu_Fsr[0] >= 110): continue
            if args.section != -1:
                if args.isdata:
                    if not (Muons_Minv_MuMu_Fsr[0] >= 110 and Muons_Minv_MuMu_Fsr[0] <= 180): continue
                    if (Muons_Minv_MuMu_Fsr[0] >= 120 and Muons_Minv_MuMu_Fsr[0] <= 130): continue
                else:
                    if not (Muons_Minv_MuMu_Fsr[0] >= 120 and Muons_Minv_MuMu_Fsr[0] <= 130): continue

                # separate events into training and validation samples
                if eventNumber[0]%4 != args.section: continue

            event_yield+=1

            # mets
            eM = []
            eM.append([metFinalTrk[0] if args.region == 'two_jet' else -9999, -9999])

            # Jets
            ej = []

            if (njet[0] == 0): 
                ej = [[-9999, -9999, -9999, -9999, -9999],[-9999, -9999, -9999, -9999, -9999]]
            elif (njet[0] == 1): 
                DeltaPhi_mumuj1[0] = ROOT.TVector2.Phi_mpi_pi(Jets_Phi_Lead[0] - Z_Phi_FSR[0])
                DeltaPhi_mumuj2[0] = -9999
                ej = [[Jets_PT_Lead[0], Jets_Eta_Lead[0], DeltaPhi_mumuj1[0], -9999, -9999],[-9999, -9999, -9999, -9999, -9999]]
            else:
                DeltaPhi_mumuj1[0] = ROOT.TVector2.Phi_mpi_pi(Jets_Phi_Lead[0] - Z_Phi_FSR[0])
                DeltaPhi_mumuj2[0] = ROOT.TVector2.Phi_mpi_pi(Jets_Phi_Sub[0] - Z_Phi_FSR[0])
                ej = [[Jets_PT_Lead[0], Jets_Eta_Lead[0], DeltaPhi_mumuj1[0], -9999, -9999],[Jets_PT_Sub[0], Jets_Eta_Sub[0], DeltaPhi_mumuj2[0], -9999, -9999]]

            # dijet
            dj = []
            if (njet[0] >= 2): 
                DeltaPhi_mumujj[0] = ROOT.TVector2.Phi_mpi_pi(Jets_Phi_jj[0]-Z_Phi_FSR[0])
                jj = ROOT.TLorentzVector()
                jj.SetPtEtaPhiM(Jets_PT_jj[0], Jets_Eta_jj[0], Jets_Phi_jj[0], Jets_Minv_jj[0])
                Jets_Y_jj[0] = jj.Rapidity()
                dj.append([Jets_PT_jj[0], Jets_Y_jj[0], DeltaPhi_mumujj[0], Jets_Minv_jj[0]])
                
            else: 
                dj.append([-9999] * 4)
                DeltaPhi_mumujj[0] = -9999

            # dimuon
            dm = []
            dm.append([Z_PT_FSR[0], Z_Y_FSR[0], -9999, abs(Muons_CosThetaStar[0])])

            met.append(eM)
            jets.append(ej)
            dijet.append(dj)
            dimuon.append(dm)

            
            # weight
            if args.section == -1:
                weight[0] = Weight_Global[0]
            else:
                wt.append(Weight_Global[0])

        if args.section == -1:
            outtree.Fill()
    
    elapsed_time = time.time() - start_time
    tt = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print '%s |%s| 100/100'%(tt, '>'*100)

    print "INFO: The event yield after preselection is: ", event_yield

    # If the arrays are not empty, then convert them into numpy arrays and save them

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
    elif 'MC16E' in args.name or 'mc16e' in args.name:
        array_basename='mc16e_'+array_basename
    elif 'data15' in args.name:
        array_basename='data15_'+array_basename
    elif 'data16' in args.name:
        array_basename='data16_'+array_basename
    elif 'data17' in args.name:
        array_basename='data17_'+array_basename
    elif 'data18' in args.name:
        array_basename='data18_'+array_basename
    else:
        print 'WARNING: File does not match with a tipical file name!!'

    region = args.region

    if event_yield==0:
        print "WARNING: No event left after the preselections. Will not save the results."
        with open('arrays/%s/%s_%s_%d.txt' % (category, array_basename, region, args.section),"w") as txt:
            txt.write('%s_%s_%d.npy\n' % (array_basename, region, args.section))
    else:
        if args.section != -1: # and not args.isdata:
            wt = np.array(wt)
            print 'Saving numpy arrays of weight to "arrays/%s"...' % (category)
            np.save('arrays/%s/%s_%s_%d_weight.npy' % (category, array_basename, region, args.section), wt)

        met = np.array(met)
        print 'Saving numpy arrays of MET to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s_%s_%d_met.npy' % (category, array_basename, region, args.section), met)

        jets = np.array(jets)
        print 'Saving numpy arrays of jets to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s_%s_%d_jets.npy' % (category, array_basename, region, args.section), jets)

        dijet = np.array(dijet)
        print 'Saving numpy arrays of dijet to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s_%s_%d_dijet.npy' % (category, array_basename, region, args.section), dijet)

        dimuon = np.array(dimuon)
        print 'Saving numpy arrays of dimuon to \"arrays/%s\"...' % (category)
        np.save('arrays/%s/%s_%s_%d_dimuon.npy' % (category, array_basename, region, args.section), dimuon)

        if args.section == -1:
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
    if not os.path.isfile("inputs/"+args.name):
        print "ERROR: The input file %s doesn\'t exist!!" % ("inputs/"+args.name)
        print 'Exiting...'
        quit()
    print "========================================================="
    print "Processing %s" % args.name
    print "========================================================="
    print ""
    print "INFO: Making %s samples" % ('general' if args.section == -1 else 'NO.' + str(args.section))

    process_arrays(args)
    
    print ""
    print "========================================================="
    print "Finishing the process"
    print "=========================================================" 
    return

if __name__ == '__main__':
    main()
