# Launch the Confinement runs
# This is program and cluster specific

import os
import pdb
import yaml
import math
# from popen2 import popen2
import subprocess
import numpy as np
from shutil import copytree, ignore_patterns, rmtree
# from sklearn.grid_search import ParameterGrid
from sklearn.model_selection import ParameterGrid


def Launch():

    # Specify options for run launcher
    options = {
            'qos' : 'condo',
            'partition' : 'shas',
            'account' : 'ucb-summit-smr',
            'time' : '72:00:00',
            # Define architecture of clusters
            'coresPerNode' : 24,
            'socksPerNode' : 2,
            }

    # Initialize the runs
    simPaths, seedPaths, nTasks = Initialize() 
    # Update Random Seeds
    UpdateRandomSeeds( ['RunConfig.yaml'], ['rngSeed'], seedPaths)

    # Write launch script
    jobStrings = WriteSeedLaunchString( options, simPaths, seedPaths)
    
    # Launch individually
    for spath, jobString in zip( seedPaths, jobStrings):
        os.chdir( spath)
        with open('amsos.jobscript', 'w') as f:
            f.write( jobString)
        # subprocess.call(["sbatch", "amsos.jobscript"])

def Initialize():
    
    # Read param file
    stream = open( "Params.yaml", 'r')
    ydict = yaml.load( stream)
    if not ydict['seeds']:
        nSeeds = 1
    else:
        nSeeds = ydict['seeds']
    
    # Get parameter grid
    grid = CreateParameterGrid( ydict)

    # Initialize directories
    simPaths, seedPaths, nTasks = InitializeDirectories( ydict, grid, nSeeds)

    # Update yaml parameters for each sim 
    UpdateYamlSim( ydict, grid, simPaths, nSeeds)

    return simPaths, seedPaths, nTasks


def CreateParameterGrid( yamldict):

    # Exit if there are no parameters to vary
    if 'parameters' not in yamldict:
        raise Exception('Parameters must be specified in Params.yaml')

    # Initialize a dictionary of parameter information
    params = {}

    # Loop over parameters and add them to the dictionary
    for key, dat in yamldict['parameters'].items():
            
        # Get Values
        values = dat["value"]

        # Get labels. Check that number of labels match es the number of values
        if len(dat["value_key"]) != len(values):
            raise Exception("length of labels for key "+str(key)+" must equal the length of values")
        labels = [ str( dat["label_key"])+str(e) for e in dat["value_key"] ]

        # Construct pairs of values and labels and add them to params
        pairs = [ (e,f) for e,f in zip(values,labels) ]
        params[key] = pairs

    # Get a parameter grid
    grid = ParameterGrid( params)
    # print('Parameters w/ labels : ')
    # for key,val in params.items():
        # print(key+" : "+str(val))
    # print('Parameter Grid : ')
    # for e in grid:
        # print(e)
    return grid

def InitializeDirectories( ydict, grid, nSeeds):
    # Initialize directories and copy all files to respective folders

    # Define parent name where the run will be stored. Initialize the dir
    if grid != ['sim']:
        parentName = ydict['dir_name']
        if os.path.exists( parentName):
            rmtree( parentName)
        os.mkdir( parentName)

    simPaths = []
    seedPaths = []
    for item in grid:
        # Construct sim path 
        if item == 'sim':
            simName = item
            simPath = os.path.join( os.getcwd(), simName) 
        else:
            simName = "_".join( [val[1] for key,val in item.items()] )
            simPath = os.path.join( os.getcwd(), parentName, simName) 
            print(simPath)

        # Delete simPath if it exists
        if os.path.exists( simPath):
            rmtree( simPath)

        # Add simPath to collection of simPaths
        simPaths += [simPath]

        for seed in range(nSeeds):

            # Construct seed folder name
            seedPath = os.path.join( simPath, 's{0}'.format(seed) )
            
            # Copy files to all seed folders
            dest = copytree( os.getcwd(), seedPath, ignore=ignore_patterns('Launch.py*',parentName, 'sim', 'data'))

            # Add seedPath to collection of seedPaths
            seedPaths += [seedPath]

    nTasks =  nSeeds * len(simPaths)
    return simPaths, seedPaths, nTasks 

def UpdateYamlSim( ydict, grid, simPaths, nSeeds):
    # Update parameters in the yaml files for each sim

    if grid == ['sim']:
        return

    print('Seeds = {0}'.format(nSeeds))
    for simPath,sim in zip( simPaths, grid):
        print('Sim : '+simPath)

        for seed in range(nSeeds):
            # Seed path
            seedPath = os.path.join( simPath, 's{0}'.format(seed) )

            # Load yaml file
            yname1 = os.path.join( seedPath,'RunConfig.yaml')
            yname2 = os.path.join( seedPath,'ProteinConfig.yaml')
            data1 = yaml.load( open( yname1, 'r') )
            data2 = yaml.load( open( yname2, 'r') )

            for key, value in sim.items():

                if key == 'boundary_diameter_tube':
                    # filament length
                    leng = data1['sylinderLength']
                    # Set tube radius 
                    # pdb.set_trace()
                    data1['boundaries'][0]['radius'] = (leng*value[0])/2

                elif key == 'packing_fraction':

                    # Get volume of system
                    height = data1['simBoxHigh'][2]- data1['simBoxLow'][2]
                    rad = data1['boundaries'][0]['radius']
                    Vol = math.pi * rad**2 * height
                    # Get volume of fil
                    leng = data1['sylinderLength']
                    rad = data1['sylinderDiameter']/2
                    vol = math.pi * rad**2 * leng

                    # num fil
                    # pdb.set_trace()
                    data1['sylinderNumber'] = int(np.floor(0.01*value[0]*Vol/vol))

            # Set protein number
            data2['freeNumber'] = int( np.ceil( ydict['protein_ratio']*data1['sylinderNumber']))

            # Write yaml file
            with open(yname1, 'w') as yaml_file:
                yaml_file.write( yaml.dump( data1))
            with open(yname2, 'w') as yaml_file:
                yaml_file.write( yaml.dump( data2))

            print('  {0} : {1}'.format(key, str(value)))

def UpdateRandomSeeds( filenames, paramnames, seedPaths):
    # Update random seed in yaml files 

    for idx, seedPath in enumerate( seedPaths):

        for fil,par in zip(filenames, paramnames):

            yname = os.path.join( seedPath,'{0}'.format(fil))
            data = yaml.load( open( yname, 'r') )

            # Update seed 
            data[par] = data[par]+idx 

            # Write yaml file
            with open(yname, 'w') as yaml_file:
                yaml_file.write( yaml.dump( data))

def WriteSeedLaunchString( slurm, simPaths, seedPaths):
    # write strings to send to sbatch for each seed. Also save it to a file in the seed directory for running later if needed

    jobStrings = []
       
    # Define jobString
    jobStringDef1 = """#!/bin/bash

#SBATCH --job-name={0}
#SBATCH --qos={1}
#SBATCH --partition={2}
#SBATCH --account={3}
#SBATCH --output=sim.log
#SBATCH --error=sim.err
#SBATCH --time={4}
#SBATCH --nodes={5}
#SBATCH --cpus-per-task={6}

export OMP_NUM_THREADS={7}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


"""
    jobStringDef2 = """#!/bin/bash

#SBATCH --job-name={0}
#SBATCH --qos={1}
#SBATCH --partition={2}
#SBATCH --account={3}
#SBATCH --output=sim.log
#SBATCH --error=sim.err
#SBATCH --time={4}
#SBATCH --nodes={5}
#SBATCH --exclusive

export OMP_NUM_THREADS={6}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


"""
    # Loop over seeds and make job strings to launch
    for spath in seedPaths:

        # Jobname : SimName_SeedNumber
        jobName = '__'.join( spath.split('/')[-2:] )
        
        # How many processors?
        yname1 = os.path.join( spath,'RunConfig.yaml')
        data1 = yaml.load( open( yname1, 'r') )
        nproc = int(np.ceil(data1['sylinderNumber']/100))

        # numThreads
        num_threads = min([nproc, slurm['coresPerNode']])
        # Write jobString 
        if nproc < slurm['coresPerNode']:
            jobString = jobStringDef1.format( jobName, slurm['qos'], slurm['partition'], 
                    slurm['account'], slurm['time'], 1, num_threads, num_threads) 
        else:
            jobString = jobStringDef2.format( jobName, slurm['qos'], slurm['partition'], 
                    slurm['account'], slurm['time'], int(np.ceil(nproc/ slurm['coresPerNode'])), num_threads) 

        jobString = jobString + 'srun -n1 --mpi=pmi2 AMSOS\n' 
        jobStrings += [jobString]

    return jobStrings

if __name__ == "__main__":
    Launch()
