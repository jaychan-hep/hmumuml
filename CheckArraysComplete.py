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
        elif 'MC16E' in filename or 'mc16e' in filename:
            #continue
            array_basename='mc16e_'+array_basename
        elif 'data15' in filename:
            array_basename='data15_'+array_basename
        elif 'data16' in filename:
            array_basename='data16_'+array_basename
        elif 'data17' in filename:
            array_basename='data17_'+array_basename
        elif 'data18' in filename:
            array_basename='data18_'+array_basename
        else:
            print 'WARNING: File does not match with a tipical file name!!'
                    
        for section in range(-1,4):
            if category in ['Z', 'stop', 'diboson', 'ttbar'] and section != -1: continue
            for region in ['two_jet','one_jet','zero_jet']:
                for obj in ['jets','dijet','dimuon','met']: 
                    if not ((os.path.isfile('arrays/%s/%s_%s_%d_%s.npy' % (category, array_basename, region, section, obj))  and (os.path.isfile('AppInputs/%s/%s_%s.root' % (category,array_basename,region)) or section != -1) and (os.path.isfile('arrays/%s/%s_%s_%d_weight.npy' % (category,array_basename, region, section)) or section == -1)) or os.path.isfile('arrays/%s/%s_%s_%d.txt' % (category,array_basename, region, section))):
#                       print '%s%s_%s_%s.npy' % (array_basename, section, region, obj)
                        samples.append([filename,str(section),region,category,channel])
                        break
    return samples

def main():
    args=getArgs()
    inputdir="inputs"
    sample_list = {'data':    ['data15','data16','data17','data18'],
                   'ttbar':   ['410472'],
                   'Z':       ['364100','364101','364102','364103','364104','364105','364106','364107','364108','364109','364110','364111','364112','364113',
                               '366300','366301','366303','366304','366306'],
                   'stop':    ['410644','410645'],
                   'diboson': ['363356','363358','364250','364253','364254'],
                   'ttH':     ['344388'],
                   'VH':      ['345103','345104','345105','345098'],
                   'VBF':     ['345106'],
                   'ggF':     ['345097']}

    samples=[]
    samples_d={}
    for category in sample_list:

        for channel in sample_list[category]:
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
        print 'Section:  %s' %('general' if sample[1]=='-1' else sample[1] )
        print 'Region:   %s' %(sample[2])
    print '====================================================================================================='
    print 'Missing summary:'
    for misch in samples_d:
        print '%s:  %d' %(misch, len(samples_d[misch]))
    print 'Total number of missing files:  %d' %(len(samples))
    print '====================================================================================================='


if __name__=='__main__':
    main()


