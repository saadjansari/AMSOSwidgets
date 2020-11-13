# Launch the runs
# This is program and cluster specific

import os, pdb, yaml
# from popen2 import popen2
import subprocess
import numpy as np
from shutil import copytree, ignore_patterns, rmtree
# from sklearn.grid_search import ParameterGrid
from sklearn.model_selection import ParameterGrid


def Launch():

    # Specify options for run launcher
    options = {
            'routine' : 'multicore', # singlecore or multicore
            'name' : 'AMSOSrun',
            'qos' : 'condo',
            'partition' : 'shas',
            'account' : 'ucb-summit-smr',
            'time' : '00:01:00',
            # Define architecture of clusters
            'coresPerNode' : 24,
            'socksPerNode' : 2,
            }

    # Initialize the runs
    simPaths, seedPaths, nTasks = Initialize() 

    # Update Random Seeds
    UpdateRandomSeeds( ['RunConfig.yaml'], ['rngSeed'], seedPaths)

    # Write launch script
    # jobStrings = WriteSimLaunchString( options, simPaths, seedPaths)
    jobStrings = WriteSeedLaunchString( options, simPaths, seedPaths)

    # # Sim Launch
    # for string in jobStrings:
        # # open pipe to sbatch command
        # output, input = popen2('sbatch')
        # # send job to sbatch
        # input.write( string)
        # input.close()
        # # Print your job and the response to the screen
        # print(string)
        # print( output.read() )
    
    # Launch individually
    for spath, jobString in zip( seedPaths, jobStrings):
        os.chdir( spath)
        with open('jobscript.sh', 'w') as f:
            f.write( jobString)
        subprocess.call(["sbatch", "jobscript.sh"])

def Initialize():
    
    # Read param file
    stream = open( "Params.yaml", 'r')
    ydict = yaml.load( stream)
    nSeeds = ydict['seeds']
    
    # Get parameter grid
    grid = CreateParameterGrid( ydict)

    # Initialize directories
    simPaths, seedPaths, nTasks = InitializeDirectories( grid, nSeeds)

    # Update yaml parameters for each sim 
    UpdateYamlSim( grid, simPaths, nSeeds)

    return simPaths, seedPaths, nTasks


def CreateParameterGrid( yamldict):

    # Exit if there are no parameters to vary
    if 'files' not in yamldict:
        return ['sim'] 

    # Initialize a dictionary of parameter information
    params = {}

    # Loop over parameters and add them to the dictionary
    for fil, fdat in yamldict['files'].items():
        for key, dat in fdat.items():
            
            # Get list of values and labels
            if dat["scanType"] == "all":

                # Get Values
                values = dat["value"]

                # Get labels. Check that number of labels match es the number of values
                if len(dat["valueKey"]) != len(values):
                    raise Exception("length of labels for key "+str(key)+" must equal the length of values")
                labels = [ str( dat["labelKey"])+str(e) for e in dat["valueKey"] ]

            elif dat["scanType"] == "interval":

                # Check integer numValues
                if type( dat["numValues"]) != int:
                    raise ValueError("numValues for key "+str(key)+" must be of type int for scanType = interval\n")

                # Get Values. Check that value is given as [lowerBound, upperBound]
                if len( dat["value"]) != 2:
                    raise ValueError("value for key "+str(key)+" must be [lowerBound, upperBound] for scanType = interval\n")
                values = np.linspace( dat["value"][0], dat["value"][1], np.ceil( dat["numValues"]) )

            else:
                raise ValueError("invalid scanType for key "+str(key)+"\n")

            # Construct pairs of values and labels and add them to params
            pairs = [ (e,f) for e,f in zip(values,labels) ]
            params[ fil+"___"+key ] = pairs

    # Get a parameter grid
    grid = ParameterGrid( params)
    # print('Parameters w/ labels : ')
    # for key,val in params.items():
        # print(key+" : "+str(val))
    # print('Parameter Grid : ')
    # for e in grid:
        # print(e)

    return grid

def InitializeDirectories( grid, nSeeds):
    # Initialize directories and copy all files to respective folders

    # Define parent name where the run will be stored. Initialize the dir
    if grid != ['sim']:
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

        # Delete simPath if it exists
        if os.path.exists( simPath):
            rmtree( simPath)

        # Add simPath to collection of simPaths
        simPaths += [simPath]

        for seed in range(nSeeds):

            # Construct seed folder name
            seedPath = os.path.join( simPath, 's{0}'.format(seed) )
            
            # Copy files to all seed folders
            dest = copytree( os.getcwd(), seedPath, ignore=ignore_patterns('Launch.py*','run*', 'sim', 'data'))

            # Add seedPath to collection of seedPaths
            seedPaths += [seedPath]

    nTasks =  nSeeds * len(simPaths)
    return simPaths, seedPaths, nTasks 

def UpdateYamlSim( grid, simPaths, nSeeds):
    # Update parameters in the yaml files for each sim

    if grid == ['sim']:
        return

    print('Seeds = {0}'.format(nSeeds))
    for simPath,sim in zip( simPaths, grid):

        print('Sim : '+simPath)
        for key, value in sim.items():

            keys = key.split('___')
            for seed in range(nSeeds):

                # Seed path
                seedPath = os.path.join( simPath, 's{0}'.format(seed) )

                # Load yaml file
                yname = os.path.join( seedPath,'{0}.yaml'.format(keys[0]))
                data = yaml.load( open( yname, 'r') )

                # Update parameter
                data[keys[1]] = value[0]

                # Write yaml file
                with open(yname, 'w') as yaml_file:
                    yaml_file.write( yaml.dump( data))

            print('  {0} : {1}'.format(key, str(value[0])))

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
    nTasks = 1
    # Find cores per socket
    coresPerSock = slurm['coresPerNode']/slurm['socksPerNode']
       
    # Define jobString
    jobStringDef = """#!/bin/bash

#SBATCH --job-name={0}
#SBATCH --qos={1}
#SBATCH --partition={2}
#SBATCH --account={3}
#SBATCH --output=sim.log
#SBATCH --error=sim.err
#SBATCH --time={4}
#SBATCH --nodes={5}
#SBATCH --cpus-per-task={6}
#SBATCH --ntasks-per-socket={7}

export OMP_NUM_THREADS={8}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


"""
    # Loop over seeds and make job strings to launch
    for spath in seedPaths:

        # Find number of nodes and number of processors/task
        if slurm['routine'] == 'singlecore':
            nCpuPerTask = 1
            nNodes = int( np.ceil(nTasks/(float(slurm['coresPerNode'])/nCpuPerTask)) )
            nTaskPerSock = coresPerSock
        elif slurm['routine'] == 'multicore':
            # raise Exception('Remove me if you want to run multicore you crazy diamond')
            nCpuPerTask = coresPerSock 
            nNodes = int( np.ceil(nTasks/(float(slurm['coresPerNode'])/nCpuPerTask)) )
            nTaskPerSock = 1

        # Jobname : SimName_SeedNumber
        jobName = '__'.join( spath.split('/')[-2:] )

        # Write jobString 
        jobString = jobStringDef.format( jobName, slurm['qos'], slurm['partition'], slurm['account'], slurm['time'], nNodes, nCpuPerTask, nTaskPerSock, nCpuPerTask) 

        jobString = jobString + 'srun -n1 --mpi=pmi2 AMSOS\n' 
        jobStrings += [jobString]

    return jobStrings

def WriteSimLaunchString( slurm, simPaths, seedPathsAll):
    # write strings to send to sbatch for each sim 

    jobStrings = []
    for spath in simPaths:

        # find seeds associated with this sim
        seedPaths = [s for s in seedPathsAll if spath in s]
    
        # Number of tasks
        nTasks = len( seedPaths)

        # Find cores per socket
        coresPerSock = slurm['coresPerNode']/slurm['socksPerNode']

        # Find number of nodes and number of processors/task
        if slurm['routine'] == 'singlecore':
            nCpuPerTask = 1
            nNodes = int( np.ceil(nTasks/(float(slurm['coresPerNode'])/nCpuPerTask)) )
            nTaskPerSock = coresPerSock
        elif slurm['routine'] == 'multicore':
            # raise Exception('Remove me if you want to run multicore you crazy diamond')
            nCpuPerTask = coresPerSock 
            nNodes = int( np.ceil(nTasks/(float(slurm['coresPerNode'])/nCpuPerTask)) )
            nTaskPerSock = 1

        # Write jobscript
        jobString = """#!/bin/bash

#SBATCH --job-name={0}
#SBATCH --qos={1}
#SBATCH --partition={2}
#SBATCH --account={3}
#SBATCH --output=sim.log
#SBATCH --error=sim.err
#SBATCH --time={4}
#SBATCH --nodes={5}
#SBATCH --cpus-per-task={6}
#SBATCH --ntasks-per-socket={7}

export OMP_NUM_THREADS={8}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


""".format( os.path.split(spath)[-1], slurm['qos'], slurm['partition'], slurm['account'], slurm['time'], nNodes, nCpuPerTask, nTaskPerSock, nCpuPerTask) 

        # For each sim_seed, add to jobString
        for pth in seedPaths:
            # command = 'srun -n1 --exclusive --mpi=pmi2 --chdir={0} AMSOS 1> {1} 2> {2} &\n'.format( pth, 'sim.log', 'sim.err')
            command = 'srun -n1 --exclusive --mpi=pmi2 --chdir={0} AMSOS 1> {1} 2> {2} &\n'.format( pth, os.path.join(pth,'sim.log'), os.path.join(pth,'sim.err'))
            jobString = jobString + command
        jobString = jobString + "wait\n"

        jobStrings += [jobString]

    return jobStrings

def WriteRunLaunchString( slurm, seedPaths):
    # write a string to send to sbatch for the entire run 

    # Number of tasks
    nTasks = len( seedPaths)

    # Find cores per socket
    coresPerSock = slurm['coresPerNode']/slurm['socksPerNode']

    # Find number of nodes and number of processors/task
    if slurm['routine'] == 'singlecore':
        nCpuPerTask = 1
        nNodes = int( np.ceil(nTasks/(float(slurm['coresPerNode'])/nCpuPerTask)) )
        nTaskPerSock = coresPerSock
    elif slurm['routine'] == 'multicore':
        raise Exception('Remove me if you want to run multicore you crazy diamond')
        nCpuPerTask = coresPerSock 
        nNodes = int( np.ceil(nTasks/(float(slurm['coresPerNode'])/nCpuPerTask)) )
        nTaskPerSock = 1

    # Write jobscript
    jobString = """#!/bin/bash

#SBATCH --job-name={0}
#SBATCH --qos={1}
#SBATCH --partition={2}
#SBATCH --account={3}
#SBATCH --output=sim.log
#SBATCH --error=sim.err
#SBATCH --time={4}
#SBATCH --nodes={5}
#SBATCH --cpus-per-task={6}
#SBATCH --ntasks-per-socket={7}

export OMP_NUM_THREADS={8}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


    """.format( slurm['name'], slurm['qos'], slurm['partition'], slurm['account'], slurm['time'], nNodes, nCpuPerTask, nTaskPerSock, nCpuPerTask) 

    # For each sim_seed, add to jobString
    for pth in seedPaths:
        # command = 'srun -n1 --exclusive --mpi=pmi2 --chdir={0} AMSOS 1> {1} 2> {2} &\n'.format( pth, 'sim.log', 'sim.err')
        command = 'srun -n1 --exclusive --mpi=pmi2 --chdir={0} AMSOS 1> {1} 2> {2} &\n'.format( pth, os.path.join(pth,'sim.log'), os.path.join(pth,'sim.err'))
        jobString = jobString + command
    jobString = jobString + "wait\n"

    return jobString


if __name__ == "__main__":
    Launch()
