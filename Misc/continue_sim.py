import os
import glob
from shutil import copytree, ignore_patterns, copyfile
import subprocess
import pdb

# Ask user for sim folder path to continue
prompt = '\nSpecify relative path to sim folder to continue: '
relpath = input(prompt) # get path from user
fpath = os.path.abspath( relpath)
if not os.path.exists( fpath):
    raise Exception('specified path does not exist')
print('Sim being continued: {}\n'.format(fpath) )

# new path contains a '_c' to signal a continuation
newpath = fpath + '_c' 

# Copy all files/folders except the 'result' folder to the newpath
dest = copytree( fpath, newpath, ignore=ignore_patterns('result')
        
# create result folder and copy few files from old resutl folder except actual results
rpath_old = os.path.join( fpath, 'result')
rpath_new = os.path.join( newpath, 'result')
dest = copytree( rpath_old, rpath_new, ignore=ignore_patterns('result*')

# copy last data files from old sim folder to use an intial files for new sim
# define a key finder
def fileKey(f):
    k = int(f[f.rfind("_")+1:f.rfind(".")])
    return k

fil1 = glob.glob( os.path.join(fpath, 'result/result*/SylinderAscii_*.dat'))
fil2 = glob.glob( os.path.join(fpath, 'result/result*/ProteinAscii_*.dat'))
fil1 = sorted(fil1, key=fileKey)
fil2 = sorted(fil2, key=fileKey)

# Copy files
shutil.copyfile( fil1[-1], os.path.join( newpath, 'TubuleInitial.dat'))
shutil.copyfile( fil2[-1], os.path.join( newpath, 'ProteinInitial.dat'))
print('Set up initial .dat files')

# Launch sim
os.chdir( newpath)
print('Launching sim...')
subprocess.call(["sbatch", "amsos.jobscript"])
