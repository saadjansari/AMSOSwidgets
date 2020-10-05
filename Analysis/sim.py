import os
import pdb
import yaml
import glob
import shutil
import numpy as np
import matplotlib.pyplot as plt
import progressbar
import pandas as pd

from decorators import timer
from Frame import *


# Class for processing a complete sim
class Sim():
    def __init__(self, folder_name, label, opts):

        self.label = label
        self.opts = opts
        self.folder_name = folder_name
        self.feather_path = os.path.join( self.folder_name, 'df.ftr') # path to a binary file storing useful processed data

    def read_config(self):
        # Read config file and store useful quantities

        # Load parameters from RunConfig.yaml and ProteinConfig.yaml files
        with open( os.path.join(self.folder_name,'RunConfig.yaml')) as f:
            pars1 = yaml.load( f, Loader=yaml.FullLoader)
        with open( os.path.join(self.folder_name,'ProteinConfig.yaml')) as f:
            pars2 = yaml.load( f, Loader=yaml.FullLoader)
            
        self.opts.dt = pars1['timeSnap'] # Specify frame step of simulation
        self.opts.boxsize = np.array( pars1['initBoxHigh']) # Specify frame step of simulation

        # Create plots directory
        if self.opts.graph:
            rpath = os.path.join(self.folder_name, 'plots')
            if os.path.exists(rpath):
                shutil.rmtree( rpath)
            os.mkdir( rpath)


    def analyze(self):
        # Analyze the sim

        # Read the config file
        self.read_config()
        print('Sim: {}'.format(self.label))

        # Load the processed binary file in case it exists
        if not self.opts.overwrite and os.path.exists( self.feather_path):
            print('\tLoading feather')
            self.df = pd.read_feather(self.feather_path)
            return 

        # If proccessed binary file does not exist, do this analysis

        self.df = pd.DataFrame()
        # Locate ascii data for filaments and proteins
        sfiles = sorted( 
                glob.glob( os.path.join(self.folder_name, 'result/result*/SylinderAscii_*.dat')),
                key=fileKey)
        pfiles = sorted( 
                glob.glob( os.path.join(self.folder_name, 'result/result*/ProteinAscii_*.dat')),
                key=fileKey)

        # Create progress bar
        bar = progress_bar( len(sfiles) )
        for idx in range( len(sfiles) ):
            
            # Initialize frame instance and launch analysis
            frame = Frame(sfiles[idx], pfiles[idx], self.opts)
            frame.analyze()
            self.df = self.df.append(frame.data, ignore_index=True)
            bar.update(idx+1)

        bar.finish()

        # save dataframe as a feather binary file
        self.df.to_feather(self.feather_path)
        print('\n')


def progress_bar(max_index):
# Create and start a progress bar
    bar = progressbar.ProgressBar(maxval=max_index, 
            widgets=[
                progressbar.Percentage(), ' ', 
                progressbar.Bar('=', '[', ']'), ' ', 
                progressbar.Timer()])
    bar.start()
    return bar

def fileKey(f):
# A function that returns the file index.
# e.g. filekey('SylinderAscii_59.dat) -> 59
# Get key to sort through file indices
    k = int(f[f.rfind("_")+1:f.rfind(".")])
    return k

if __name__ == "__main__":

    fpath = "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin_6400/run/f1"

    sim = Sim( fpath, 'f1', opts)
    sim.analyze()
    pdb.set_trace()
    print('woah buddy')

