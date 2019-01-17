#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  8.13.2018
#
#
#
#
#

import os
from datetime import datetime
from argparse import ArgumentParser

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for MonoHbb XGBoost analysis.")
    parser.add_argument('-i', '--inputdir', action='store', default='inputs', help='Directory that contains ntuples.')
    return  parser.parse_args()

def submitLSF(params,jobname,date):

    initDir=os.getcwd()
    lsfDir=initDir+"/condor/SubmitProcessArrays/"+date+"/"+jobname+"/"
    if not os.path.exists(lsfDir):
       os.makedirs(lsfDir)

    jdl="#Agent jdl file\n"
    jdl+="Universe        = vanilla\n"
    jdl+="Notification    = Never\n"
    jdl+="initialdir      = "+initDir+"\n"
    jdl+="Executable      = "+initDir+"/SubmitRecoverMissingFiles.sh\n"
    #jdl+="GetEnv          = True\n"
    jdl+="Output          = "+lsfDir+jobname+".$(ClusterId).$(ProcId).out\n"
    jdl+="Error           = "+lsfDir+jobname+".$(ClusterId).$(ProcId).err\n"
    jdl+="Log             = "+lsfDir+jobname+".$(ClusterId).$(ProcId).log\n"
    jdl+="stream_output   = False\n"
    jdl+="stream_error    = False\n"
    #jdl+="should_transfer_files = yes\n"
    #jdl+='Requirements = ((Arch == "X86_64") && (regexp("CentOS",OpSysLongName)))\n'
    jdl+='Requirements = ((Arch == "X86_64") && (regexp("SLC",OpSysLongName)))\n'
    jdl+="WhenToTransferOutput = ON_EXIT_OR_EVICT\n"
    jdl+="OnExitRemove         = TRUE\n"
    jdl+='+JobFlavour = "tomorrow"\n'
    jdl+='+JobType="HyperTuning"\n'
    jdl+='+AccountingGroup ="group_u_ATLASWISC.all"\n'
    jdl+="RequestCpus = 1\n"
    #jdl+="Arguments = "+params+" "+initDir+" \nQueue \n"
    for param in params:
        jdl+="Arguments = "+initDir+" "+param+" \nQueue \n"

    jdlFile=lsfDir+jobname+".jdl"
    handle=open(jdlFile,"w")
    handle.write(jdl)
    handle.close()
    command="chmod +x "+jdlFile
    os.system(command)
    if jdlFile==None:
       print "JDL is None\n"
       sys.exit(1)

    command="condor_submit "+jdlFile
    print command
    os.system(command)
    return




def main():
    args=getArgs()
    date = datetime.now().strftime("%Y-%m-%d-%H-%M")
    inputdir=args.inputdir
    sample_list=['data','ttH','VH','VBF','ggF','Z2']

    samples={}
    for channel in sample_list:
        #if channel=='Znunu': continue
        samples[channel]=[]
    
    #print samples

    text_file = open("arrays/missing_samplelist.txt", "r")
    lines = text_file.read().splitlines()
    for line in lines:
        args=str(line).split(' ')
        if 'train' in args[1]:
            section='t'
        elif 'val' in args[1]:
            section='v'
        else:
            section='none'
        #print samples
        #print line
        samples[args[3]]+=["%s %s %s %s"% (args[0],section,args[2],args[3])]
        #print samples
        #break

        
    for channel in samples:
        if not samples[channel]==[]: submitLSF(samples[channel],channel,date)

                


    print 'Log files will be saved to \"condor/SubmitProcessArrays/%s/\"'%(date)
    print 'Use \"condor_q\" to check to status of condor jobs.'
    
if __name__=='__main__':
    main()
