import os
import pdb
import glob
import subprocess

path1 = '/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf'
scripts = glob.glob( os.path.join(path1, '*/*/*/Result2PVD.py'))

for spath in scripts:

    pth = '/'.join( spath.split('/')[:-1])
    print(pth)
    os.chdir( pth)
    subprocess.call(["python", "Result2PVD.py"])
