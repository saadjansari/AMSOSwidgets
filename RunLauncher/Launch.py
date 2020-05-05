# Launch the runs
# This is program and cluster specific

import os, pdb, yaml
import subprocess
import numpy as np
from shutil import copytree, ignore_patterns, rmtree
# from sklearn.grid_search import ParameterGrid
from sklearn.model_selection import ParameterGrid

def Initialize():
    
    # Read param file
    stream = open( "Params.yaml", 'r')
    ydict = yaml.load( stream)
    
    # Get parameter grid
    grid = CreateParameterGrid( ydict)

    # Initialize directories
    folds = InitializeDirectories( grid)

    # Update Parameters
    UpdateParameters( grid, folds)

def UpdateParameters( grid, folds):
    # Update parameters in the run yaml files 

    for fold,run in zip( folds, grid):

        print('Run :'+fold)
        for key, value in run.items():

            keys = key.split('___')
            # Load yaml file
            yname = os.path.join( fold, keys[0] + '.yaml')
            data = yaml.load( open( yname, 'r') )

            # Update parameter
            data[keys[1]] = value[0]

            # Write yaml file
            with open(yname, 'w') as yaml_file:
                # yaml_file.write( yaml.dump( data, default_flow_style=False))
                yaml_file.write( yaml.dump( data))

            print('  '+key+' : '+str(value[0]))

def InitializeDirectories( grid):
    # Initialize directories and copy all files to respective folders

    parentName = 'runs'
    if os.path.exists( parentName):
        rmtree( parentName)
    os.mkdir( parentName)

    folds = []
    for item in grid:

        # Construct folder name
        labs = [val[1] for key,val in item.items()]
        folderName = "_".join(labs)
        fpath = os.path.join( os.getcwd(), parentName, folderName)

        # Copy all files to folder
        if os.path.exists( fpath):
            rmtree( fpath)
        dest = copytree( os.getcwd(), fpath, ignore=ignore_patterns('Launch.py*','Initialize.py*', 'run*'))
        folds += [fpath]

    return folds

def CreateParameterGrid( yamldict):

    # Initialize a dictionary of parameter information
    params = {}

    # Loop over parameters and add them to the dictionary
    for fil, fdat in yamldict.items():
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

def WriteLaunchScript():
    # create a bash script for executing the jobs

    nJobs = WriteJobsFile()
    
    f_launch = open('run.launch.sh', "w+")
    # identiy bash evironment and define sbatch options
    f_launch.write( '#!/bin/bash\n\n')
    f_launch.write( '#SBATCH --job-name=TubSim\n')
    f_launch.write( '#SBATCH --qos=condo\n')
    f_launch.write( '#SBATCH --partition=shas\n')
    f_launch.write( '#SBATCH --account=ucb-summit-smr\n')
    f_launch.write( '#SBATCH --output=sim.log\n')
    f_launch.write( '#SBATCH --error=sim.err\n')
    f_launch.write( '#SBATCH --time=00:05:00\n')

    routine = 'singlecore'
    coresPerNode = 24
    socksPerNode = 2
    coresPerSock = 12
    if routine == 'singlecore':
        nCpuPerTask = 1
        nNodes = int( np.ceil(nJobs/(float(coresPerNode)/nCpuPerTask)) )
        f_launch.write( '#SBATCH --nodes={0}\n'.format(nNodes) )
        f_launch.write( '#SBATCH --ntasks={0}\n'.format(nJobs) )
        f_launch.write( '#SBATCH --cpus-per-task={0}\n'.format(nCpuPerTask) )

    elif routine == 'multicore':
        raise Exception('Remove me if you want to run multicore you crazy diamond')
        nCpuPerTask = socksPerNode 
        nNodes = int( np.ceil(nJobs/(float(coresPerNode)/nCpuPerTask)) )
        f_launch.write( '#SBATCH --nodes={0}\n'.format(nNodes) )
        f_launch.write( '#SBATCH --ntasks={0}\n'.format(nJobs) )
        f_launch.write( '#SBATCH --ntasks-per-socket=1\n')
        f_launch.write( '#SBATCH --cpus-per-task={0}\n'.format(nCpuPerTask) )

    f_launch.write( '\nmodule load loadbalance\n')

    f_launch.write( '\nexport OMP_NUM_THREADS={0}\n'.format(nCpuPerTask))
    f_launch.write( 'export OMP_PROC_BIND=spread\n')
    f_launch.write( 'export OMP_PLACES=threads\n')

    # Job execution command
    f_launch.write( '\nmpirun lb {0}\n'.format('run.jobs') )

    # close the file
    f_launch.close()
    print( 'Successful write to launch file!') 

def WriteJobsFile():
    # write list of bash commands to be executed by loadbalancer

    # open file
    print( 'Writing jobs file...')
    f_jobs = open('run.jobs', "w+")

    # Prefix for job command
    cmd_pre = 'srun --mpi=pmi2 '
    # Suffix for job command
    cmd_post = '"\n'
    
    # Get list of runs to start
    folds = os.listdir( os.path.join(os.getcwd(), 'runs'))
    
    for fpath in folds:
        cmd_launch = 'cd {0}; srun --mpi=pmi2 AMSOS \n'.format( os.path.join(os.getcwd(),'runs',fpath) )
        f_jobs.write( cmd_launch) 
    
    f_jobs.close()
    print( 'Successful write to jobs file!') 

    return len(folds) 

if __name__ == "__main__":

    # Initialize the runs
    Initialize() 

    # Launch individually
    os.chdir( 'runs')
    cwd = os.getcwd()
    files = os.listdir( os.getcwd())
    for fil in files:
        os.chdir( os.path.join( cwd,fil))
        subprocess.call(["sbatch", "amsos.jobscript"])

    # Write launch script
    # WriteLaunchScript()
