#!/usr/bin/env python

import sys, os, pdb
import glob, re 
import pickle, copy, shutil
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from unit_dict import UnitDict
import Plots

ud = UnitDict()

class Sim(object):

    def __init__(self, cwd=None, label='default'):

        if not cwd:
            self.cwd = os.getcwd() 
        else:
            self.cwd = cwd

        self.seeds = []
        self.pickleFile = 'simData.pickle'
        self.label = label

        self.InitSeeds()


    def InitSeeds(self):

        # Specify path of pickle file
        fname = os.path.join( self.cwd, self.pickleFile)

        # If pickle file exists, load it
        if os.path.exists( fname):
            print('Loading sim pickle file...')
            with open( fname, 'rb') as f:

                # Load it
                self.frames, self.config = pickle.load(f)

            # Update the working directory of loaded seeds 
            for seed in self.seeds:
                seed.cwd = os.path.join( self.cwd, os.path.split( seed.cwd)[-1])
                # Create seed directories if they dont exist
                if not os.path.exists( seed.cwd) :
                    os.mkdir( seed.cwd)

        # If pickle file doesn't exist, initialize seeds
        else:

            # Find seed directories
            seedList = glob.glob( self.cwd+'/s*')

            # Initialize each seed 
            for f in seedList:
                print('   Seed: '+ f)
                self.seeds += [Seed( os.path.join( self.cwd,f), f)]

            # Get config 
            self.config = self.seeds[0].config

            # Save pickle file
            print('Saving sim pickle file...')
            with open(fname, 'wb') as f:
                pickle.dump( [self.frames, self.config], f) 


    def GraphQuantity(self, quantity, ax=None):

        # Create figure, ax if not provided
        if not ax:
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )

        # Get timestep 
        timestep = self.config['RunConfig']['timeSnap']

        # Plot for each seed
        cols = sns.color_palette('husl', len(self.seeds))
        for seed,c in zip(self.seeds, cols):
            seed.GraphQuantity( quantity, ax=ax, color=c, label=seed.label):

    def GraphQuantityMean(self, quantity, ax=None):

        # Create figure, ax if not provided
        if not ax:
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )

        # Get timestep 
        timestep = self.config['RunConfig']['timeSnap']

        # Plot for each seed
        cols = sns.color_palette('husl', len(self.seeds))
        for seed,c in zip(self.seeds, cols):
            seed.GraphQuantity( quantity, ax=ax, color=c, label=seed.label):

    def Graph(self, quantities=None):

        # Initialize graphing data directory
        gdir = os.path.join( self.cwd, 'data')
        if os.path.exists( gdir):
            shutil.rmtree( gdir, ignore_errors=True)

        # Specify quantities to graph
        if not quantities:
            quantities = ['polarOrder, nematicOrder']

        os.chdir( gdir)
        for quant in quantities:

            # Create figure 
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )
            self.GraphQuantity( quant, ax)
            fig.suptitle( quant)
            fig.savefig( quant+'.pdf')
            plt.close()
        os.chdir( self.cwd)


if __name__ == '__main__':
    print('not implemented yet')    
