#!/usr/bin/python

#module needed python

import os
from glob import glob

for file in glob("./*merged*.fa"):
    file = file.strip()
    (x, name) = file.split('/', 1)
    (name, x) = name.split('.', 1)
    fa2kerJob = "fa2ker_{0}".format(name)
    fa2kerOut = 'fa2ker_{0}.out'.format(name)
    Output = "{0}_get2ker.csv".format(name)
    Command = "sbatch -t 1-00:00:00 --mem=80000 -J {0} -o {1} --wrap=\"python fa2kmerAbundanceTable.py -i {2} -o {3} -k 2 --percentage\"".format(fa2kerJob, fa2kerOut, file, Output)
    #print(Command)
    os.system(Command)
