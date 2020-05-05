#!/usr/bin/env python

import sys, os, pdb
import glob, re 
import pickle, yaml, vtk, copy, shutil
import argparse
# from sklearn.grid_search import ParameterGrid
from sklearn.model_selection import ParameterGrid
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from unit_dict import UnitDict
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

    opts = parser.parse_args()
    return opts

# member variables are dynamically added by parsing data files
# for Sylinder
class Sylinder(object):
    end0 = None
    end1 = None
    pass

class Frame(object):
    def __init__(self, sylinderFile=None):
        self.sylinders = []
        self.parseSylinderFile(sylinderFile)

    def parseFile(self, dataFile, objType, objList):
        # print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
        for i in range(nObj):
            s = objType()
            s.end0 = data.GetPoints().GetPoint(2 * i)
            s.end1 = data.GetPoints().GetPoint(2 * i + 1)
            objList.append(s)

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            # print("Parsing Cell Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName, cdata.GetTuple(j))

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            # print("Parsing Point Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName + "0", pdata.GetTuple(2 * j))
                setattr(objList[j], dataName + "1", pdata.GetTuple(2 * j + 1))

    def parseSylinderFile(self, sylinderFile):
        self.parseFile(sylinderFile, Sylinder, self.sylinders)


class Sim(object):

    def __init__(self, cwd=None, label='default'):

        if not cwd:
            self.cwd = os.getcwd() 
        else:
            self.cwd = cwd

        self.frames = []
        self.pickleFile = 'simData.pickle'
        self.label = label

        # If pickle file exists, load it
        pklFull = os.path.join( self.cwd, self.pickleFile)
        if os.path.exists( pklFull):
            print('Loading sim pickle file...')
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
            print('Saving sim pickle file...')
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

    def Graph(self):
        print('1')

    def PlotSylinder(self):
        print('1')

class Run(object):

    def __init__(self, cwd = None, label = 'default'):

        if not cwd:
            self.cwd = os.getcwd() 
        else:
            self.cwd = cwd

        self.sims= []
        self.params = []
        self.pickleFile = 'runData.pickle'
        self.label = label

        # If pickle file exists, load it
        pklFull = os.path.join( self.cwd, self.pickleFile)
        if os.path.exists( pklFull):
            print('Loading run pickle file...')
            with open( pklFull, 'rb') as f:
                self.sims, self.params = pickle.load(f)
            # Create sim directories if they dont exist
            for sim in self.sims:
                if not os.path.exists( sim.cwd) :
                    os.mkdir( sim.cwd)
        else:
            # get list of sim directories. filter directories starting with a '.'
            simFileList = next(os.walk(self.cwd))[1]
            simFileList = [i for i in simFileList if not i[0] == '.']
            # initialize each sim
            for fil in simFileList:
                print('   Sim: '+ fil)
                self.sims += [Sim( os.path.join( self.cwd,fil), fil)]

            # Load params
            with open( os.path.join(self.sims[0].cwd,'Params.yaml')) as f:
                self.params = yaml.load( f)

            # Save pickle file
            print('Saving run pickle file...')
            with open(pklFull, 'wb') as f:
                pickle.dump( [self.sims,self.params], f) 

        # Specify variables to graph
        self.graph_vars_onedog = ['velBrown', 'velNonBrown']


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
        

    def Graph(self):
        
        if os.path.exists(  os.path.join( self.cwd, 'data')):
            shutil.rmtree( os.path.join( self.cwd,'data'), ignore_errors=True)
        os.mkdir( os.path.join( self.cwd,'data'))
        os.chdir( os.path.join(self.cwd,'data'))

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

            fig.suptitle(figname)
            fig.savefig(figname+'.pdf')
            plt.close()


class Analysis(object):
    def __init__(self, opts):

        self.opts = opts
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



def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


if __name__ == '__main__':
    
    opts = parseArgs()
    ana = Analysis( opts)

