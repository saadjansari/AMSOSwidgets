import os
import pdb
import glob
import shutil

path1 = os.path.join( os.getcwd(), 'eq')
path1 = os.path.join( os.getcwd(), 'run')

def fileKey(f):
    k = int(f[f.rfind("_")+1:f.rfind(".")])
    return k

# eq/run folder has specific sim directories, which have seed directories
# look in those seed directories
sfil = glob.glob( os.path.join(path1, '*/*'))
for spath in sfil:

    fil1 = glob.glob( os.path.join(spath, 'result/result*/SylinderAscii_*.dat'))
    fil2 = glob.glob( os.path.join(spath, 'result/result*/ProteinAscii_*.dat'))
    fil1 = sorted(fil1, key=fileKey)
    fil2 = sorted(fil2, key=fileKey)

    # path to copy to
    cpath = os.path.join( path2, '/'.join(spath.split('/')[-2:]) )

    # Copy files
    print('Copied from {0} to {1}'.format(fil1[-1], os.path.join( cpath, 'TubuleInitial.dat') ))
    shutil.copyfile( fil1[-1], os.path.join( cpath, 'TubuleInitial.dat'))
    shutil.copyfile( fil2[-1], os.path.join( cpath, 'ProteinInitial.dat'))
