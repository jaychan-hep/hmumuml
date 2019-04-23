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

def submitLSF(params,jobname,date):

    initDir=os.getcwd()
    lsfDir=initDir+"/condor/cat_opt/"+date+"/"
    if not os.path.exists(lsfDir):
       os.makedirs(lsfDir)

    jdl="#Agent jdl file\n"
    jdl+="Universe        = vanilla\n"
    jdl+="Notification    = Never\n"
    jdl+="initialdir      = "+initDir+"\n"
    jdl+="Executable      = "+initDir+"/submit_categorization_optimization_2D.sh\n"
    #jdl+="GetEnv          = True\n"
    jdl+="Output          = "+lsfDir+jobname+".$(ClusterId).$(ProcId).out\n"
    jdl+="Error           = "+lsfDir+jobname+".$(ClusterId).$(ProcId).err\n"
    jdl+="Log             = "+lsfDir+jobname+".$(ClusterId).$(ProcId).log\n"
    jdl+="stream_output   = False\n"
    jdl+="stream_error    = False\n"
    #jdl+="should_transfer_files = yes\n"
    #jdl+='Requirements = ((Arch == "X86_64") && (regexp("CentOS",OpSysLongName)))\n'
    jdl+='Requirements = ((Arch == "X86_64") && (regexp("CentOS7",OpSysAndVer)))\n'
    #jdl+='Requirements = ((Arch == "X86_64") && (regexp("SLC",OpSysLongName)))\n'
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
    date = datetime.now().strftime("%Y-%m-%d-%H-%M")

    arg=[]
    for i in range(2,8):
        for j in range(2,8):
            #if i!=j: continue
            #for k in range(135,165):
            for k in range(50,90):
                if os.path.isfile('cat_opt/%d_%d/two_jet_0_%d.txt'%(j,i,k)): continue
                arg.append('%d %d %d'%(k, i, j))
    submitLSF(arg,'cat_opt',date)

    print 'Log files will be saved to \"condor/cat_opt/%s/\"'%(date)
    print 'Use \"condor_q\" to check to status of condor jobs.'

if __name__=='__main__':
    main()

