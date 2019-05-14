#!/usr/bin/env python
#
#
#
#  Created by Jay Chan
#
#  3.13.2019
#
#
#
#
#
import os
from datetime import datetime
from argparse import ArgumentParser
from ttHyy.condor import condor_booklist

def main():

    date = datetime.now().strftime("%Y-%m-%d-%H-%M")

    condor_list = condor_booklist('submit_categorization_optimization_2D.sh', 'cat_opt')
    condor_list.initialdir_in_arguments()
    condor_list.set_JobFlavour('longlunch')

    for i in range(2,7):
        for j in range(2,7):
            if i!=j: continue
            #for k in range(80,180):
            for k in range(40,90):
                if os.path.isfile('cat_opt/two_jet/%d_%d/%d.json'%(j,i,k)): continue
                condor_list.add_Argument('%d %d %d'%(k, i, j))

    condor_list.submit()

if __name__=='__main__':
    main()

