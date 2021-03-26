#!/usr/bin/env python
#
#
#  Created by Jay Chan (jay.chan@cern.ch)
#
#     05.05.2019
#
#
#
import os
from datetime import datetime

def createScript(path, script):

    dir_path = os.path.dirname(path)

    if dir_path and not os.path.isdir(dir_path):
        os.makedirs(dir_path)

    file = open(path, "w")
    file.write(script)
    file.close()

    return

class condor_booklist(object):
    """Class for condor job bookkeeping"""

    def __init__(self, Executable = '', JobType = 'default_condor', JobName = ''):
        self.Universe = 'vanilla'
        self.Notification = 'Never'
        self.initialdir = os.getcwd()
        self.Executable = Executable
        self.GetEnv = ''
        self.stream_output = False
        self.stream_error = False
        self.should_transfer_files = ''
        self.Requirements = '((Arch == "X86_64") && (regexp("CentOS7",OpSysAndVer)))'
        self.WhenToTransferOutput = 'ON_EXIT_OR_EVICT'
        self.OnExitRemove = True
        self.JobFlavour = 'tomorrow'
        self.JobType = JobType
        self.AccountingGroup = 'group_u_ATLASWISC.all'
        self.RequestCpus = 1
        self.Request_memory = "2000 MB"
        self.Request_disk = "1 KB"
        self.Arguments = []
        self.initialdir_argument = False
        self.JobName = JobName

    def set_Universe(self, Universe):
        self.Universe = Universe

    def set_Notification(self, Notification):
        self.Notification = Notification

    def set_initialdir(self, initialdir):
        self.initialdir = initialdir

    def set_Executable(self, Executable):
        self.Executable = Executable

    def set_GetEnv(self, GetEnv):
        self.GetEnv = GetEnv

    def set_stream_output(self, stream_output):
        self.stream_output = stream_output

    def set_stream_error(self, stream_error):
        self.stream_error = stream_error

    def set_should_transfer_files(self, should_transfer_files):
        self.should_transfer_files = should_transfer_files

    def set_Requirements(self, Requirements):
        self.Requirements = Requirements

    def set_WhenToTransferOutput(self, WhenToTransferOutput):
        self.WhenToTransferOutput = WhenToTransferOutput

    def set_OnExitRemove(self, OnExitRemove):
        self.OnExitRemove = OnExitRemove

    def set_JobFlavour(self, JobFlavour):
        self.JobFlavour = JobFlavour

    def set_JobType(self, JobType):
        self.JobType = JobType

    def set_AccountingGroup(self, AccountingGroup):
        self.AccountingGroup = AccountingGroup

    def set_RequestCpus(self, RequestCpus):
        self.RequestCpus = RequestCpus

    def set_Request_memory(self, Request_memory):
        self.Request_memory = Request_memory

    def set_Request_disk(self, Request_disk):
        self.Request_disk = Request_disk

    def set_JobName(self, JobName):
        self.JobName = JobName

    def set_Arguments(self, Arguments):
        self.Arguments = Argements

    def add_Argument(self, Argument):
        if isinstance(Argument, list):
            self.Arguments.extend(Argument)
        else:
            self.Arguments.append(Argument)

    def add_BlankArgument(self):
        self.add_Argument(' ')

    def initialdir_in_arguments(self):
        self.initialdir_argument = True

    def summary(self,keyword=''):
        print('----------Condor Booklist Sumary----------')
        if not keyword or keyword == 'Basic' or keyword == 'Executable':
            print('Executable:          {Executable}'.format(Executable = self.Executable))
        if not keyword or keyword == 'Basic' or keyword == 'Executable' or keyword == 'initialdir':
            print('Initial directory:   {initialdir}'.format(initialdir = self.initialdir))
        if not keyword or keyword == 'Basic' or keyword == 'JobType':
            print('Job type:            {JobType}'.format(JobType = self.JobType))
        if not keyword or keyword == 'Basic' or keyword == 'JobType' or keyword == 'JobName' and self.JobName:
            print('Job name:            {JobName}'.format(JobName = self.JobName))
        if not keyword or keyword == 'Basic' or keyword == 'JobFlavour':
            print('Job Favour:          {JobFlavour}'.format(JobFlavour = self.JobFlavour))
        if not keyword or keyword == 'Basic' or keyword == 'Arguments' and self.Arguments:
            print('-----Arguments-----')
            for argument in self.Arguments:
                print(argument)
            print('----------------')
        print('------------------------------------------')

    def submit(self, with_arguments = True):

        if not self.Executable:
            print('ERROR: executable is not specified!!')
            quit()

        if not self.Arguments:
            if not with_arguments: Self.Arguments.Append(' ')
            else:
                print('WARNING: no jobs are added. Will not submit any job.')
                return

        date = datetime.now().strftime("%Y-%m-%d-%H-%M")        
        lsfDir = '{initialdir}/condor/{JobType}/{date}{JobName}'.format(initialdir = self.initialdir, JobType = self.JobType, date = date, JobName = '' if not self.JobName else '.' + self.JobName)
        if not os.path.exists(lsfDir):
            os.makedirs(lsfDir)

        jdl =  "#Agent jdl file\n"
        jdl += "Universe        = {Universe}\n".format(Universe = self.Universe)
        jdl += "Notification    = {Notification}\n".format(Notification = self.Notification)
        jdl += "initialdir      = {initialdir}\n".format(initialdir = self.initialdir)
        jdl += "Executable      = {initialdir}/{Executable}\n".format(initialdir = self.initialdir, Executable = self.Executable)
        if self.GetEnv: jdl += "GetEnv          = {GetEnv}\n".format(GetEnv = self.GetEnv)
        jdl += "Output          = {lsfDir}/$(ClusterId).$(ProcId).out\n".format(lsfDir = lsfDir)
        jdl += "Error           = {lsfDir}/$(ClusterId).$(ProcId).err\n".format(lsfDir = lsfDir)
        jdl += "Log             = {lsfDir}/$(ClusterId).$(ProcId).log\n".format(lsfDir = lsfDir)
        jdl += "stream_output   = {stream_output}\n".format(stream_output = self.stream_output)
        jdl += "stream_error    = {stream_error}\n".format(stream_error = self.stream_error)
        if self.should_transfer_files: jdl += "should_transfer_files = {should_transfer_files}\n".format(should_transfer_files = self.should_transfer_files)
        jdl += 'Requirements    = {Requirements}\n'.format(Requirements = self.Requirements)
        jdl += "WhenToTransferOutput = {WhenToTransferOutput}\n".format(WhenToTransferOutput = self.WhenToTransferOutput)
        jdl += "OnExitRemove         = {OnExitRemove}\n".format(OnExitRemove = self.OnExitRemove)
        jdl += '+JobFlavour = "{JobFlavour}"\n'.format(JobFlavour = self.JobFlavour)
        jdl += '+JobType="{JobType}"\n'.format(JobType = self.JobType)
        if self.AccountingGroup: jdl += '+AccountingGroup ="{AccountingGroup}"\n'.format(AccountingGroup = self.AccountingGroup)
        jdl += "RequestCpus = {RequestCpus}\n".format(RequestCpus = self.RequestCpus)
        jdl += "Request_memory = {Request_memory}\n".format(Request_memory = self.Request_memory)
        jdl += "Request_disk = {Request_disk}\n".format(Request_disk = self.Request_disk)
        for Argument in self.Arguments:
            jdl += "Arguments = {initialdir} {Argument} \nQueue \n".format(initialdir = self.initialdir if self.initialdir_argument else '', Argument = Argument)

        jdlFile = "{lsfDir}/{JobName}.jdl".format(lsfDir = lsfDir, JobName = self.JobName if self.JobName else 'default_condor')
        handle = open(jdlFile, "w")
        handle.write(jdl)
        handle.close()
        command = "chmod +x " + jdlFile
        os.system(command)
        if jdlFile == None:
            print("JDL is None\n")
            sys.exit(1)

        command = "condor_submit " + jdlFile
        print(command)
        os.system(command)

        return
