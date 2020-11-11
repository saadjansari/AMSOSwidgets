# Launch the Confinement runs
# This is program and cluster specific

import os
import pdb
import yaml
import math
import subprocess
import glob
import numpy as np
import argparse
from shutil import copytree, ignore_patterns, rmtree, copyfile
# from sklearn.grid_search import ParameterGrid
from sklearn.model_selection import ParameterGrid


'''
Name: Launch.py
Description: Construct and launch a run with multiple sims with specifically constrcuted parameter files
Input: To see type Launch.py -h
'''

def parseArgs():

    parser = argparse.ArgumentParser(prog='Launch.py')
    parser.add_argument('--sbatch', action='store_true', default=False, help='Submit jobs via bash sbatch to cluster.')
    parser.add_argument('-E', '--eq', action='store_true', default=False, help='Launch equilibration')
    parser.add_argument('-I', '--use_eq_init', action='store_true', default=False, help='Use initial files present in a similar directory tree.')

    opts = parser.parse_args()

    return opts


def Launch(opts):

    # Initialize the runs
    simPaths, seedPaths, nTasks = Initialize(opts) 

    # Write launch script
    jobStrings = WriteSeedLaunchString( simPaths, seedPaths)
    
    # Launch individually
    for spath, jobString in zip( seedPaths, jobStrings):
        os.chdir( spath)
        with open('amsos.jobscript', 'w') as f:
            f.write( jobString)
        if opts.sbatch:
            subprocess.call(["sbatch", "amsos.jobscript"])
            

def Initialize(opts):
    
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
    simPaths, seedPaths, nTasks = InitializeDirectories( ydict, grid, nSeeds,opts)

    # Update yaml parameters for each sim 
    UpdateYamlSim( ydict, grid, simPaths, nSeeds, opts)

    # Update Random Seeds
    UpdateRandomSeeds( ['RunConfig.yaml'], ['rngSeed'], seedPaths)

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
    return grid

def InitializeDirectories( ydict, grid, nSeeds, opts):
    # Initialize directories and copy all files to respective folders

    # Define parent name where the run will be stored. Initialize the dir
    if grid != ['sim']:
        if 'dir_name' in ydict.keys():
            parentName = ydict['dir_name']
        elif opts.eq:
            parentName = 'eq'
        else:
            parentName = 'run'

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
            dest = copytree( os.getcwd(), seedPath, ignore=ignore_patterns('Launch.py*',parentName, '*eq', 'eq*', '*run', 'run*', 'sim', 'data'))

            # Add seedPath to collection of seedPaths
            seedPaths += [seedPath]

    nTasks =  nSeeds * len(simPaths)

    # Copy init equilibrium files if specified
    if opts.use_eq_init:
        # eq path
        p1 = os.path.join( os.getcwd(), 
                str( input("Specify relative path of parent directory of equilibrium run : ") ) )
        # Get run path
        p2 = os.path.join( os.getcwd(), parentName)

        # eq/run folder has specific sim directories, which have seed directories
        # look in those seed directories
        sfil = glob.glob( os.path.join(p1, '*/*'))
        for spath in sfil:
            fil1 = glob.glob( os.path.join(spath, 'result/result*/SylinderAscii_*.dat'))
            fil2 = glob.glob( os.path.join(spath, 'result/result*/ProteinAscii_*.dat'))
            fil1 = sorted(fil1, key=fileKey)
            fil2 = sorted(fil2, key=fileKey)

            # path to copy to
            cpath = os.path.join( p2, '/'.join(spath.split('/')[-2:]) )

            # Copy files
            print('Copied from {0} to {1}'.format(fil1[-1], os.path.join( cpath, 'TubuleInitial.dat') ))
            shutil.copyfile( fil1[-1], os.path.join( cpath, 'TubuleInitial.dat'))
            shutil.copyfile( fil2[-1], os.path.join( cpath, 'ProteinInitial.dat'))
        
    return simPaths, seedPaths, nTasks 

def UpdateYamlSim( ydict, grid, simPaths, nSeeds, opts):
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

            if ydict['regime'] == 'Tactoid':
                data1, data2 = create_tactoid_config( sim, ydict, data1,data2)
                if opts.eq:
                    data1['sylinderFixed'] = True
            elif ydict['regime'] == 'Confinement':
                data1, data2 = create_confinement_config( sim, ydict, data1,data2)
                if opts.eq:
                    data2['proteins'][0]['Ka'] = [0,0]

            # Write yaml file
            with open(yname1, 'w') as yaml_file:
                yaml_file.write( yaml.dump( data1))
            with open(yname2, 'w') as yaml_file:
                yaml_file.write( yaml.dump( data2))


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
                

def WriteSeedLaunchString( simPaths, seedPaths):
    # write strings to send to sbatch for each seed. Also save it to a file in the seed directory for running later if needed

    # ask user for time
    val_time = str( input( "Specify job time HH:MM:SS : ") )

    # Specify options for run launcher
    slurm = {
            'qos' : 'condo',
            'partition' : 'shas',
            'account' : 'ucb-summit-smr',
            'time' : val_time,
            # Define architecture of clusters
            'coresPerNode' : 24,
            'socksPerNode' : 2,
            }

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
    # Loop over seeds and make job strings to launch
    for spath in seedPaths:

        # Jobname : SimName_SeedNumber
        jobName = '__'.join( spath.split('/')[-2:] )
        
        # How many processors?
        yname1 = os.path.join( spath,'RunConfig.yaml')
        data1 = yaml.load( open( yname1, 'r') )
        nproc = min( [int(np.ceil(data1['sylinderNumber']*0.01)) , slurm['coresPerNode'] ])

        # numThreads
        num_threads = min([nproc, slurm['coresPerNode']])
        # Write jobString 
        jobString = jobStringDef1.format( jobName, slurm['qos'], slurm['partition'], 
                slurm['account'], slurm['time'], 1, num_threads, num_threads) 

        jobString = jobString + 'srun -n1 --mpi=pmi2 AMSOS\n' 
        jobStrings += [jobString]

    return jobStrings

def create_confinement_config( sim, ydict, data1, data2):

    # Parameters in Param.yaml and what they vary
    # packing_fraction: varies the packing fraction by setting sylinderNumber, the number of filaments. It does this by first finding the Volume of the system and a single filament, and then calculating the number of filaments

    # boundary_diameter_tube: varies the radius of the "tube" boundary. When radius of the confining boundary is changed, the system volume changes and therefore the packing fraction automatically becomes higher. We allow a variable cylinder height to keep packing fraction constant

    # vol of filament
    leng = data1['sylinderLength']
    rad = data1['sylinderDiameter']/2
    vol = math.pi * rad**2 * leng

    # Set tube radius 
    if 'boundary_diameter_tube' in sim.keys():
        tube_rad = sim['boundary_diameter_tube'][0]*leng/2
        data1['boundaries'][0]['radius'] = tube_rad

    # packing fraction
    if 'packing_fraction' in sim.keys():
        pf = sim['packing_fraction'][0]

        if ydict['parameters']['boundary_diameter_tube']['vary_height'] is True:

            # Get volume of system
            heights = ydict['parameters']['boundary_diameter_tube']['heights']
            values = ydict['parameters']['boundary_diameter_tube']['value']
            height = heights[np.argwhere(sim['boundary_diameter_tube'][0] == np.array(values))[0][0]]*leng
            data1['simBoxHigh'][2] = data1['simBoxLow'][2]+height
            data1['initBoxHigh'][2] = data1['initBoxLow'][2]+height

            rad = data1['boundaries'][0]['radius']
            Vol = math.pi * tube_rad**2 * height

            data1['sylinderNumber'] = int(np.floor(0.01*pf*Vol/vol))

        else:
            # Get volume of system
            height = data1['simBoxHigh'][2]- data1['simBoxLow'][2]
            rad = data1['boundaries'][0]['radius']
            Vol = math.pi * tube_rad**2 * height
            # Get volume of fil
            leng = data1['sylinderLength']
            rad = data1['sylinderDiameter']/2
            vol = math.pi * rad**2 * leng

            data1['sylinderNumber'] = int(np.floor(0.01*pf*Vol/vol))

    # Set protein number
    data2['proteins'][0]['freeNumber'] = int( np.ceil( ydict['conf_params']['protein_ratio']*data1['sylinderNumber']))

    return data1,data2

def create_tactoid_config( sim, ydict, data_rc, data_pc):
    # data_rc is data from RunConfig.yaml file
    # data_pc is data from ProteinConfig file

    # Parameters in Param.yaml and what they vary
    # filamin_conc: concentration of filamin xlinker in percentage of actin concentration (10% represents 7 xlinkers per 1 actin filament)

    if 'filamin_conc' in sim.keys():
        fc = sim['filamin_conc'][0]

        # num xlinkers
        n_actin = data_rc['sylinderNumber']
        n_xlink = fc*0.7*n_actin

        # calculate how many fixed xlinks and how many free
        num_fixed = int( n_xlink // n_actin )
        num_free = int( n_xlink % n_actin )
        
        data_pc['proteins'][0]['freeNumber'] = num_free
        data_pc['proteins'][0]['fixedLocationPerMT'] = [2 for ii in range(num_fixed)]

    if 'filamin_ratio' in sim.keys():
        fc = int( sim['filamin_ratio'][0] )
        num_fixed = int( fc )
        data_pc['proteins'][0]['fixedLocationPerMT'] = [2 for ii in range(num_fixed)]

    if 'filamin_len' in sim.keys():
        fl = sim['filamin_len'][0]
        data_pc['proteins'][0]['freeLength'] = fl
        
    if 'filamin_Ka' in sim.keys():
        ka = sim['filamin_Ka'][0]
        data_pc['proteins'][0]['Ka'] = ka 

    return data_rc,data_pc

def fileKey(f):
    k = int(f[f.rfind("_")+1:f.rfind(".")])
    return k

if __name__ == "__main__":
    opts = parseArgs()
    Launch(opts)
