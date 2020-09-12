#!/usr/bin/env python3

import sys, os, pdb
import glob, re 
import pickle, yaml, vtk, copy, shutil
import argparse
from sklearn.model_selection import ParameterGrid
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from unit_dict import UnitDict
from AMSOSrun import Run 
from AMSOSsim import Sim
from AMSOSseed import Seed 
import Plots

ud = UnitDict()

'''
Name: AnalyzeAMSOS.py
Description: Command that calls specific bulk analysis programs
Input: To see type AnalyzeAMSOS.py -h
Output: See above
'''

def parseArgs():

    parser = argparse.ArgumentParser(prog='AnalyzeAMSOS.py')
    parser.add_argument('-R', '--run', action='store_true', default=False, help='Launch analysis of runs')
    parser.add_argument('-S', '--sim', action='store_true', default=False, help='Launch analysis of sims')
    parser.add_argument('-G', '--graph', action='store_true', default=False, help='Make graphs')
    parser.add_argument('--clean', action='store_true', default=False, help='Clean pickle files and exit')


    opts = parser.parse_args()
    return opts

class Analysis(object):
    def __init__(self, opts):

        self.opts = opts

        # Clear
        if opts.clean:
            print('Cleaning pickle files')
            # XXX : Python 3
            # files = glob.iglob('**.pickle', recursive=True)
            # pdb.set_trace()
            # for f in files:
                # os.remove( f)
            # XXX : Python 2
            files = glob.glob('./*.pickle') + glob.glob('./*/*.pickle') + glob.glob('./*/*/*.pickle')
            for f in files:
                os.remove( f)
            return

        if opts.run:
            print('Analyzing Run')
            run = Run()
            if opts.graph:
                run.Graph()

        elif opts.sim:
            print('Analyzing Sim')
            sim = Sim()
            if opts.graph:
                sim.Graph()


if __name__ == '__main__':
    
    opts = parseArgs()
    ana = Analysis( opts)
