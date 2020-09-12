#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pdb
import math

import os
import glob

outputName = 'opcf'

def fileKey(f):
    k = int(f[f.rfind("_")+1:f.rfind(".")])
    return k


class Tubule:
    def __init__(self, linestring):
        data = linestring.split()
        self.gid = int(data[1])
        self.minus = np.array([float(data[2]), float(data[3]), float(data[4])])
        self.plus = np.array([float(data[5]), float(data[6]), float(data[7])])
        xi = (self.plus-self.minus)
        self.orientation = xi/np.sqrt(xi.dot(xi))


class Cell:
    def __init__(self, linestring):
        data = linestring.split()
        self.gid = int(data[1])
        self.radius = float(data[2])
        self.minus = np.array([float(data[3]), float(data[4]), float(data[5])])
        self.plus = np.array([float(data[6]), float(data[7]), float(data[8])])
        xi = (self.plus-self.minus)
        self.orientation = xi/np.sqrt(xi.dot(xi))

class Frame:
    def __init__(self, filename):
        self.TList = []
        self.PList = []
        file = open(filename, 'r')
        for line in file:
            if line.startswith('T '):  # parse tubule data
                self.TList.append(Tubule(linestring=line))
            elif line.startswith('C '):  # parse 'cell' data, also put into TList
                self.TList.append(Cell(linestring=line))
            elif line.startswith('P '):
                self.PList.append(Protein(linestring=line))

def calcRDF(frame):
    # calc radial dist function 

    rmax = 100*frame.TList[0].radius
    rr = np.linspace( 0, rmax, num=60)
    dr = 4*frame.TList[0].radius
    gr = []
    opcf1 = []
    opcf2 = []
    for rc in rr:
        grs = []
        op1s = []
        op2s = []
        for tub in frame.TList:
            p0 = (tub.plus+tub.minus)/2
            for tubx in frame.TList:
                px = (tubx.plus+tubx.minus)/2
                if np.absolute( np.linalg.norm( p0-px) -rc ) <= dr:
                    grs.append( [1])
                    op1s.append( [np.dot( tub.orientation, tubx.orientation)/( np.linalg.norm(tub.orientation) * np.linalg.norm(tubx.orientation) ) ])
                    op2s.append( [0.5*( (3*np.dot( tub.orientation, tubx.orientation)/( np.linalg.norm(tub.orientation) * np.linalg.norm(tubx.orientation) ))**2 -1)] )
                else:
                    grs.append( [0])
                    op1s.append( [0])
                    op2s.append( [0])
        gr.append( np.mean(grs))
        opcf1.append( np.mean(op1s))
        opcf2.append( np.mean(op2s))

    fig, ax = plt.subplots( 1, 1)
    ax.plot( rr, gr)
    fig.savefig('rdf.pdf')

    fig, ax = plt.subplots( 1, 1)
    ax.plot( rr, opcf1)
    fig.savefig('opcf1.pdf')

    fig, ax = plt.subplots( 1, 1)
    ax.plot( rr, opcf2)
    fig.savefig('opcf2.pdf')



def opcf(frame):
    # calc orientation
    orientList = []
    for T in frame.TList:
        orientList.append(T.orientation)
    # mean
    print('Entries in frame: ', len(frame.TList))
    PList = np.array(orientList)
    QList = np.array([np.outer(p, p) for p in PList])
    polarOrder = np.mean(PList, axis=0)
    nematicOrder = np.mean(QList, axis=0) - np.identity(3)/3
    # This is the correct S
    S = np.sqrt(np.tensordot(nematicOrder, nematicOrder)*1.5)
    return np.hstack([S, polarOrder, nematicOrder.flatten()])


fil = 'SylinderAscii_999.dat'
# files = glob.glob('SylinderAscii_*.dat')
# files = sorted(files, key=fileKey)
# print('Parsing '+str(len(files))+' SylinderAscii files')

data = []

frame = Frame(fil)
rdf = calcRDF( frame)
    # print('S: \n', S)
    # print('polar order: \n', polarOrder)
    # print('nematic order: \n', nematicOrder)
    # data.append(frame_data)

# data = np.array(data)
# np.savetxt('record_'+outputName+'.csv', data, fmt='%6g', delimiter=',',
           # header='S,px,py,pz,Qxx,Qxy,Qxz,Qyx,Qyy,Qyz,Qzx,Qzy,Qzz')

# print("mean S: \n", np.mean(data[-meanWindow:, 0]))
# print("mean polar order: \n", np.mean(data[-meanWindow:, 1:4], axis=0))
# print("mean nematic order: \n", np.mean(data[-meanWindow:, 4:], axis=0))

# plt.plot(np.linalg.norm(data[:, 1:4], axis=1), 'r-', label=r'$|p|$')
# plt.plot(data[:, 0], 'b-',
         # label=r'$S=<P_2(\cos\theta)>=\sqrt{\frac{3}{2}Q:Q}$')
# plt.legend()
# plt.savefig('record_'+outputName+'.png')
# plt.show()
