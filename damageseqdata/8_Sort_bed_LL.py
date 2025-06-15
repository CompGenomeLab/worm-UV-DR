#Sort bed files
#module needed python

import os
from glob import glob

for file in glob("./*_10.bed"):
    file = file.strip()
    (x, name) = file.split('/', 1)
    (name, x) = name.split('.', 1)
    Bed_Job = 'sort_{0}'.format(name)
    Bed_Out = 'sort_{0}.out'.format(name)
    Output = '{0}_sort.bed'.format(name)
    Command = "sbatch -t 1-00:00:00 --mem=20000 -J {0} -o {1} --wrap=\"sort -k1,1 -k2,2n {2} > {3}\"".format(Bed_Job, Bed_Out, file, Output)
    os.system(Command)
