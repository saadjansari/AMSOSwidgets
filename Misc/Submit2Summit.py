import os
import pdb
import glob
import subprocess

path1 = os.path.join( os.getcwd(), 'run')
scripts = glob.glob( os.path.join(path1, '*/*/amsos.jobscript'))

for spath in scripts:

    os.chdir( '/'.join( spath.split('/')[:-1]) )
    subprocess.call(["sbatch", "amsos.jobscript"])
