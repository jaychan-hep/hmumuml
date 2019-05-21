import os
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-m','--recovermissing', action='store_true', default=False, help='Recover the missing files.')
parser.add_argument('-b','--recoverbad', action='store_true', default=False, help='Recover the bad files.')
args=parser.parse_args()

def recover(recoverlist):
    text_file = open(recoverlist, "r")
    lines = text_file.read().splitlines()
    for line in lines:
        arg=str(line).split(' ')
        if arg[1] == '-1':
            section=-1
        else:
            section=arg[1]
        print "python process_arrays.py -n %s -s %s -r %s -c %s %s" % (arg[0],section,arg[2],arg[3], '-d' if 'data' in arg[3] else '')
        os.system("python process_arrays.py -n %s -s %s -r %s -c %s %s" % (arg[0],section,arg[2],arg[3], '-d' if 'data' in arg[3] else ''))


if args.recovermissing:
    print "Recovering missing files..."
    recover('arrays/missing_samplelist.txt')
    print "All missing samples have been recovered. Will delete the missing samplelist."
    os.remove('arrays/missing_samplelist.txt')
if args.recoverbad:
    print "Recovering bad files..."
    recover('arrays/badsamplelist.txt')
    print "All missing samples have been recovered. Will delete the bad samplelist."
    os.remove('arrays/badsamplelist.txt')
if not (args.recovermissing or args.recoverbad):
    print "Recovering missing files..."
    recover('arrays/missing_samplelist.txt')
    print "All missing samples have been recovered. Will delete the missing samplelist."
    os.remove('arrays/missing_samplelist.txt')
    print "Recovering bad files..."
    recover('arrays/badsamplelist.txt')
    print "All missing samples have been recovered. Will delete the bad samplelist."
    os.remove('arrays/badsamplelist.txt')


