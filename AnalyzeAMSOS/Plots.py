#!/usr/bin/python

import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
from unit_dict import UnitDict
ud = UnitDict()

def main(frames):

    fig, axs = plt.subplots( 1,2,figsize=(12,6))

    plotSylinderAllQuantityNorm( frames, 'velBrown', axs[0])
    plotSylinderMeanQuantityNorm( frames, 'velBrown', axs[1])
    plt.show()

def plotSylinderAllQuantityNormTime( frames, quant, ax, timestep=1):
    # input: frames is a list of objects of class Frame
    
    for i in range(len(frames[0].sylinders)):
        vec = [np.linalg.norm( getattr( frame.sylinders[i], quant) ) for frame in frames]
        ax.plot( vec, timestep*np.array( range(len(frames))) )

    ax.set_xlabel( 'time (s)')
    ax.set_ylabel( quant)

def plotSylinderMeanQuantityNormTime( frames, quant, ax, col, label='', timestep=1):
    # input: frames is a list of objects of class Frame
   
    quantArr = np.zeros( (len(frames[0].sylinders), len(frames)))
    for i in range(len(frames[0].sylinders)):
        quantArr[i,:] = [np.linalg.norm( getattr( frame.sylinders[i], quant) ) for frame in frames]

    times = timestep*np.array( range(len(frames)))
    means = np.mean( quantArr, axis=0)
    stds = np.std( quantArr, axis=0)

    ax.plot(times, means, color=col)
    ax.fill_between(times, means-stds, means+stds, color=col, alpha=0.3)
    ax.set_xlabel( 'time '+ud.ud['time'])
    ax.set_ylabel( quant + ' '+ud.ud[quant])


if __name__ == "__main__":
    main()
