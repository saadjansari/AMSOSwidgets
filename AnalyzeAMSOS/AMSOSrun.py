#!/usr/bin/env python

import sys, os, pdb
import glob, re 
import pickle, vtk, copy, shutil
from sklearn.model_selection import ParameterGrid
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from unit_dict import UnitDict
from AMSOSsim import Sim
import Plots

ud = UnitDict()

class Run(object):

    def __init__(self, cwd = None, label = 'default'):

        if not cwd:
            self.cwd = os.path.join( os.getcwd() )
        else:
            self.cwd = cwd

        self.sims= []
        self.params = []
        self.pickleFile = 'runData.pickle'
        self.label = label

        self.InitSims()

        # Specify variables to graph
        # self.graph_vars_onedog = ['velBrown']
        self.graph_vars_onedog = ['velBrown', 'force', 'nematicOrder']


    def InitSims(self):

        # Specify path of pickle file
        fname = os.path.join( self.cwd, self.pickleFile)

        # If pickle file exists, load it
        if os.path.exists( fname):
            print('Loading run pickle file...')
            with open( fname, 'rb') as f:

                # Load it
                self.sims, self.params = pickle.load(f)

            # Update the working directory of loaded sims and seeds
            for sim in self.sims:
                sim.cwd = os.path.join( self.cwd, os.path.split( sim.cwd)[-1])
                # Create sim directories if they dont exist
                if not os.path.exists( sim.cwd) :
                    os.mkdir( sim.cwd)
                # Create seed directories if they dont exist
                for seed in sim.seeds:
                    seed.cwd = os.path.join( sim.cwd, os.path.split( seed.cwd)[-1])
                    if not os.path.exists( seed.cwd):
                        os.mkdir( seed.cwd)

        # If pickle file doesn't exist, initialize sims 
        else:

            # Find sim directories
            simList = next(os.walk(self.cwd))[1]
            simList = [i for i in simList if not i[0] == '.']

            # Initialize each sim 
            for f in simList:
                print('   Sim: '+ f)
                self.sims += [Sim( os.path.join( self.cwd,f), f)]

            # Load params
            with open( os.path.join(self.sims[0].cwd,'Params.yaml')) as f:
                self.params = yaml.load( f)

            # Save pickle file
            print('Saving run pickle file...')
            with open(fname, 'wb') as f:
                pickle.dump( [self.sims,self.params], f) 


    def ProcessParams(self):

        varabs = []
        labels = []
        for key, value in self.params.items():
            varabs += [ (key,i) for i in value.keys()]
            labels += [i['labelKey'] for i in self.params['RunConfig'].values()]

        values = {}
        for key in varabs:
            values[key[0]+'_'+key[1]] = {}
        for sim in self.sims:
            for key in varabs:
                values[key[0]+'_'+key[1]][sim.label] = sim.config[key[0]][key[1]]

        # Get unique values of each key
        uniques = {}
        print( 'Variable parameters:')
        for key in varabs:
            uniques[key[0]+'_'+key[1]] = np.unique( np.array( [values[key[0]+'_'+key[1]].values()]) )
            print( '  {0} : {1}'.format(key[0]+'_'+key[1], uniques[key[0]+'_'+key[1]]) )

        # Get all folders corresponding to unique param values
        uniqueFolds = {}
        for skey in varabs:
            jkey = '_'.join(skey)
            uniqueFolds[jkey] = []
            for i in range(len(uniques[jkey])):
                uniqueFolds[jkey] += [ [sim.label for sim in self.sims if sim.config[skey[0]][skey[1]] == uniques[ jkey][i]] ]

        return varabs, uniques, uniqueFolds 
        

    def Graph(self, quantities=['polarOrder, nematicOrder']):
        
        # Initialize graphing data directory
        gdir = os.path.join( self.cwd, 'data')
        if os.path.exists( gdir):
            shutil.rmtree( gdir, ignore_errors=True)
        os.chdir( gdir)

        # Process parameters for graphing
        varabs, uq, uqf = self.ProcessParams()

        # For each free parameter, get a grid for all other parameters
        for keys in varabs:
            keyj = '_'.join(keys)
            varabsX = copy.copy( varabs)
            varabsX.remove( keys)
            
            # Remove the free parameter and only keep the parameters to fix
            uqfc = copy.deepcopy( uqf)
            uqfc.pop(keyj)

            # Make a grid for all combinations of the fixed parameters
            grid = ParameterGrid( uqfc)
            for entry in grid:

                # Get an intersection of all the element sets of this entry 
                folds = set.intersection( *[ set( v) for v in entry.values()] )

                # Get sims corresponding to this list of labels
                sims = [ sim for sim in self.sims if any( [sim.label == fold for fold in folds]) ]

                # Call graphing function
                self.GraphParamOneDog( sims, free_var=keys, fixed_var=varabsX)

        # Graph Sims
        for sim in self.sims:
            sim.Graph( quantities)

            os.chdir( self.cwd)
    
    def GraphParamOneDog( self, sims, free_var, fixed_var):

        # Create directories
        if not os.path.exists(  os.path.join( self.cwd, 'data')):
            os.mkdir( os.path.join(self.cwd,'data'))
        os.chdir( os.path.join(self.cwd,'data'))

        # graphing for each variable
        for var in self.graph_vars_onedog:
            if not os.path.exists(  os.path.join( self.cwd, 'data', var)):
                os.mkdir( os.path.join( self.cwd,'data', var))
            os.chdir( os.path.join(self.cwd,'data', var))
            
            # Create figure 
            fig, axs = plt.subplots( 1,1, figsize=(6,6) )
            figname = '__'.join( [ vary[0]+'_'+vary[1]+'_'+ str(sims[0].config[vary[0]][vary[1]]) for vary in fixed_var])

            clrs = sns.color_palette('husl', len(sims))
            for sim, c in zip(sims, clrs):

                # Get timestep and label
                timestep = sim.config['RunConfig']['timeSnap']
                lab = free_var[1]+' = '+str(sim.config[free_var[0]][free_var[1]])
                Plots.plotSylinderMeanQuantityNormTime( sim.frames, var, axs, c, label=lab, timestep=timestep )
                
            if isinstance( axs, np.ndarray):
                for ax in axs:
                    ax.legend()
            else:
                axs.legend()

            fig.suptitle(figname)
            fig.savefig(figname+'.pdf')
            plt.close()
