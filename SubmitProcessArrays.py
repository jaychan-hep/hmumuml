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
    parser = ArgumentParser(description="Process input rootfiles into numpy arrays for Hmumu XGBoost analysis.")
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
    jdl+="Executable      = "+initDir+"/SubmitProcessArrays.sh\n"
    #jdl+="GetEnv          = True\n"
    jdl+="Output          = "+lsfDir+jobname+".$(ClusterId).$(ProcId).out\n"
    jdl+="Error           = "+lsfDir+jobname+".$(ClusterId).$(ProcId).err\n"
    jdl+="Log             = "+lsfDir+jobname+".$(ClusterId).$(ProcId).log\n"
    jdl+="stream_output   = False\n"
    jdl+="stream_error    = False\n"
    #jdl+="should_transfer_files = yes\n"
    #jdl+='Requirements = ((Arch == "X86_64") && (regexp("CentOS",OpSysLongName)))\n'
    jdl+='Requirements = ((Arch == "X86_64") && (regexp("CentOS7",OpSysAndVer)))\n'
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




def arrange_samplelist(channel,category,inputdir):
    samples=[]
    for filename in os.listdir(inputdir):
        if channel not in filename: continue
        if '.root' not in filename: continue
        #if 'mc16e' in filename: continue
        for section in ['g', '0', '1', '2', '3']:
            for region in ['two_jet','one_jet','zero_jet']:
                samples.append("%s %s %s %s"% (filename,category,section,region))
    return samples

def main():
    args=getArgs()
    date = datetime.now().strftime("%Y-%m-%d-%H-%M")
    inputdir=args.inputdir
    sample_list=[['data',['data15','data16','data17','data18']],
                 #['ttbar',['410472']],
                 #['Z2',['361107',
                 #      '308093']],
                 #['stop',['410015','410016']],
                 #['diboson',['363356','363358','364250','364253','364254']],
                 ['ttH',['344388']],
                 ['VH',['345103','345104','345105','345098']],
                 ['VBF',['345106']],
                 ['ggF',['345097']]]

    for j in range(len(sample_list)):
        category=sample_list[j][0]
        #if not (category=='stop' or category=='VH' or category=='diboson'): continue
        for channel in sample_list[j][1]:
            #if channel=='Znunu': continue
            samples=arrange_samplelist(channel,category,inputdir)
            submitLSF(samples,channel,date)
            #break
        #break
    
    print 'Log files will be saved to \"condor/SubmitProcessArrays/%s/\"'%(date)
    print 'Use \"condor_q\" to check to status of condor jobs.'
    
if __name__=='__main__':
    main()

