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
from argparse import ArgumentParser

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('--checkonly', action='store_true', default=False, help='Not save the text file')
    return  parser.parse_args()


def missing_samplelist(channel,category,inputdir):
    samples=[]
    for filename in os.listdir(inputdir):
        if channel not in filename:
            continue
        if not ".root" in filename:
            continue
        DSID=filename.split('.')[1]
        array_basename=DSID
        if 'MC16A' in filename or 'mc16a' in filename:
            array_basename='mc16a_'+array_basename
        elif 'MC16D' in filename or 'mc16d' in filename:
            array_basename='mc16d_'+array_basename
        elif 'data15' in filename:
            array_basename='data15_'+array_basename
        elif 'data16' in filename:
            array_basename='data16_'+array_basename
        elif 'data17' in filename:
            array_basename='data17_'+array_basename
        else:
            print 'WARNING: File does not match with a tipical file name!!'

                    
        for section in ['_train','_val','']:
            for region in ['zero_jet','one_jet','two_jet']:
                for obj in ['met','jets','muons','dijet','dimuon','extras']: 
                    if not ((os.path.isfile('arrays/%s/%s%s_%s_%s.npy' % (category,array_basename, section, region, obj))  and (os.path.isfile('AppInputs/%s/%s_%s.root' % (category,array_basename,region)) or section != '') and (os.path.isfile('arrays/%s/%s%s_%s_weight.npy' % (category,array_basename, section, region)) or section=='')) or os.path.isfile('arrays/%s/%s%s_%s.txt' % (category,array_basename, section, region))):
#                       print '%s%s_%s_%s.npy' % (array_basename, section, region, obj)
                        samples.append([filename,section,region,category,channel])
                        break
    return samples

def main():
    args=getArgs()
    inputdir="inputs"
    sample_list=[['data',['data15','data16','data17']],
                 #['ttbar',['410472']],
                 #['Z2',['361107',
                 #      '308093']],
                 #['stop',['410015','410016']],
                 #['diboson',['363356','363358','364250','364253','364254']],
                 ['ttH',['344388']],
                 ['VH',['345103','345104','345105','345098']],
                 ['VBF',['345106']],
                 ['ggF',['345097']]]

    samples=[]
    samples_d={}
    for j in range(len(sample_list)):
        category=sample_list[j][0]
        for channel in sample_list[j][1]:
            print 'Checking %s...' %(channel)
            mislist=missing_samplelist(channel,category,inputdir)
            if mislist==[]: continue
            samples+=mislist
            samples_d[channel]=mislist


    # save to txt
    if not args.checkonly:
        outfile = open('arrays/missing_samplelist.txt', 'w')
        for sample in samples:
            item=''
            for lab in sample:
                item+=(lab+' ')
            outfile.write("%s\n" % item)

    print '====================================================================================================='
    print ' Missing samples:'
    print '====================================================================================================='
    j=0
    for sample in samples:
        if j: 
            print '-----------------------------------------------------------------------------------------------------'
        else:
            j=1
        print 'File:     %s' %(sample[0])
        print 'Section:  %s' %('general' if sample[1]=='' else sample[1].replace('_','') )
        print 'Region:   %s' %(sample[2])
    print '====================================================================================================='
    print 'Missing summary:'
    for misch in samples_d:
        print '%s:  %d' %(misch, len(samples_d[misch]))
    print 'Total number of missing files:  %d' %(len(samples))
    print '====================================================================================================='


if __name__=='__main__':
    main()


