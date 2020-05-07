#!/usr/bin/env python

import sys, os, pdb
import glob, re 
import pickle, copy, shutil
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from unit_dict import UnitDict
from AMSOSseed import Seed 
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
                self.seeds, self.config = pickle.load(f)

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
                pickle.dump( [self.seeds, self.config], f) 

    def Graph(self, quantities=['polarOrder, nematicOrder']):

        # Initialize graphing data directory
        gdir = os.path.join( self.cwd, 'data')
        if os.path.exists( gdir):
            shutil.rmtree( gdir, ignore_errors=True)

        os.chdir( gdir)
        for quant in quantities:

            # Create figure 
            fig, ax = plt.subplots( 1,2, figsize=(12,6) )
            self.GraphSeedsProperty( quant, ax[0])
            self.GraphSeedsPropertyMean( quant, ax[1])

            fig.suptitle( '{0}__{1}'.format(sim.label,quant))
            fig.savefig( quant+'.pdf')
            plt.close()

        # Graph Seeds
        for seed in self.seeds:
            seed.Graph( quantities)

        os.chdir( self.cwd)


    def GraphSeedsProperty(self, quantity, ax=None, color=None, label=None):

        # Create figure, ax if not provided
        if not ax:
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )
        if not color:
            color = sns.color_palette('husl', len(self.seeds))
        if not label:
            label = [s.label for s in self.seeds]

        # Get timestep 
        timestep = self.config['RunConfig']['timeSnap']

        # Plot for each seed
        for seed,c,lab in zip(self.seeds, color, label):
            seed.GraphProperty( quantity, ax=ax, color=c, label=seed.label)
        ax.legend()

    def GraphSeedsPropertyMean(self, quantity, ax=None, color='m', label=None, err_alpha=0.3):

        # Create figure, ax if not provided
        if not ax:
            fig, ax = plt.subplots( 1,1, figsize=(6,6) )

        # Get timestep 
        timestep = self.config['RunConfig']['timeSnap']

        # Data for each seed
        means = np.zeros( (len(self.seeds), len(self.seed[0].frames)) )
        stds = np.zeros( (len(self.seeds), len(self.seed[0].frames)) )
        for i, seed in enumerate( self.seeds):
            means[i,:], std[i,:] = Plots.getSylinderMeanPropertyNorm( seed.frames, quantity)

        # Find mean and std
        mean = np.mean( means, axis=0)
        std = np.mean( stds, axis=0)
        time = timestep*np.arange( len(mean) )

        # Plot 
        ax.plot(time, mean, color=color)
        ax.fill_between(time, mean-std, mean+std, color='m', alpha=err_alpha)
        ax.set_xlabel( Plots.getPropertyLabel('time') )
        ax.set_ylabel( Plots.getPropertyLabel(quantity) )


if __name__ == '__main__':
    print('not implemented yet')    
