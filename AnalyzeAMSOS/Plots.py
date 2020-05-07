#!/usr/bin/python

import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
from unit_dict import UnitDict
ud = UnitDict()

def getSylinderPropertyNorm( frames, quant):
    # input: frames is a list of objects of class Frame
    
    nSyl = len( frames[0].sylinders)
    propAll = np.zeros( (nSyl, len(frames)))
    for i in range( nSyl):
        propAll[i,:] = [np.linalg.norm( getattr( frame.sylinders[i], quant) ) for frame in frames]
    return propAll 

def getSylinderMeanPropertyNorm( frames, quant):
    # input: frames is a list of objects of class Frame
    
    propAll = getSylinderPropertyNorm( frames, quant)
    shape = propAll.shape()

    times = timestep*np.arange( len(shape[0]) )
    means = np.mean( quantArr, axis=0)
    stds = np.std( quantArr, axis=0)

    return means, stds

def getPropertyLabel( quant):
    # Get label with units
    if not ud.ud[quant]:
        lab = quant
    else:
        lab = '{0} {1}'.format(quant, ud.ud[quant])
    return lab

def graphSylinderAllPropertyNorm( frames, quant, ax, timestep=1, alpha=1):
    # input: frames is a list of objects of class Frame
    
    propAll = getSylinderPropertyNorm( frames, quant)
    for row in propAll:
        ax.plot( row, timestep*np.arange(len(row)), alpha=alpha)

    ax.set_xlabel( getPropertyLabel('time') )
    ax.set_ylabel( getPropertyLabel(quant) )

def plotSylinderMeanPropertyNorm( frames, quant, ax, color='m', label='', timestep=1, err_alpha=0.3):
    # input: frames is a list of objects of class Frame
   
    means, stds = getSylinderMeanPropertyNorm( frames, quant)
    times = timestep*np.arange( len(means) )

    ax.plot(times, means, color=color, label=label)
    ax.fill_between(times, means-stds, means+stds, color=col, alpha=err_alpha)
    ax.set_xlabel( getPropertyLabel('time') )
    ax.set_ylabel( getPropertyLabel(quant) )


if __name__ == "__main__":
    main()
