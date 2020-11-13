#!/usr/bin/env python3
from pathlib import Path
from os import chdir
from subprocess import call

# This script launches the local 'Result2PVD.py' script for AMSOS sims. This can be used with a run path

# Prompt the user for the path to the run directory
prompt = '\nSpecify path to run folder with sims/seeds: '
relpath = Path(input(prompt))  # get path from user
fpath = Path(relpath).resolve()
if not fpath.exists():
    raise Exception('specified path does not exist')
print('Run path: {}\n'.format(fpath))
print('Creating PVD files for Paraview...')

rpvd_paths = fpath.glob('**/Result2PVD.py')
for spath in rpvd_scripts:
    spath.resolve()
    chdir( spath.parent)
    call(["python", "Result2PVD.py"])
