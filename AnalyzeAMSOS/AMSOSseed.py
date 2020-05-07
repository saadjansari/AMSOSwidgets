#!/usr/bin/env python

import sys, os, pdb
import glob, re 
import pickle, yaml, copy, shutil
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from unit_dict import UnitDict
import Plots
from ReadVTK import *

ud = UnitDict()

class Seed(object):

    def __init__(self, cwd=None, label='default'):

        if not cwd:
            self.cwd = os.getcwd() 
        else:
            self.cwd = cwd

        self.frames = []
        self.pickleFile = 'seedData.pickle'
        self.label = label

        # If pickle file exists, load it
        pklFull = os.path.join( self.cwd, self.pickleFile)
        if os.path.exists( pklFull):
            print('Loading seed pickle file...')
            with open( pklFull, 'rb') as f:
                self.frames, self.config = pickle.load(f)
        else:
            # get file list, sort as numerical order
            SylinderFileList = glob.glob( self.cwd+'/result/result*/Sylinder_*.pvtp')
            SylinderFileList.sort(key=getFrameNumber_lambda)
            # Get frames
            for i in range(len(SylinderFileList)):
                self.frames += [Frame(SylinderFileList[i])]
            # Load config
            self.LoadConfig()
            # Save pickle file
            print('Saving seed pickle file...')
            with open(pklFull, 'wb') as f:
                pickle.dump( [self.frames, self.config], f) 

    def LoadConfig( self):
        
        self.config = {
                'RunConfig': [],
                'ProteinConfig':[]
                }
        # load runconfig and proteinconfig yaml files
        with open(self.cwd+'/RunConfig.yaml') as f:
            self.config['RunConfig'] = yaml.load( f)
        with open(self.cwd+'/ProteinConfig.yaml') as f:
            self.config['ProteinConfig'] = yaml.load( f)

    def GraphProperty(self, quantity, ax=None, color='m', label=None):

        # Create figure, ax if not provided
        if not ax:
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )

        # Get timestep 
        timestep = self.config['RunConfig']['timeSnap']
        Plots.graphSylinderMeanPropertyNorm( self.frames, quantity, ax, color=color, timestep=timestep)

    def Graph(self, quantities=['polarOrder, nematicOrder']):

        # Initialize graphing data directory
        gdir = os.path.join( self.cwd, 'data')
        if os.path.exists( gdir):
            shutil.rmtree( gdir, ignore_errors=True)
        os.chdir( gdir)

        for quant in quantities:

            # Create figure 
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )
            self.GraphProperty( quant, ax)
            fig.suptitle( '{0}__{1}'.format(seed.label,quant))
            fig.savefig( quant+'.pdf')
            plt.close()

        os.chdir( self.cwd)

def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))

if __name__ == '__main__':
    print('not implemented yet')    
