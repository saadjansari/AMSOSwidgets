#!/usr/bin/env python3

import sys, os, pdb
import glob, re 
import vtk 
import numpy as np

class Sylinder(object):
    end0 = None
    end1 = None
    orientation = None
    pass

class Frame(object):
    def __init__(self, sylinderFile=None):
        self.sylinders = []
        self.parseSylinderFile(sylinderFile)
        self.calcOrderParameters()

    def calcOrderParameters(self):
        # Calculate nematic and order parameters for sylinders

        # Get orientations
        orients = []
        for sylinder in self.sylinders:
            orients.append( sylinder.orientation)
        orients = np.array( orients)
         
        # Polar order
        polarOrder = np.mean( orients, axis=0)

        # Nematic order tensors
        Qs = np.array( [np.outer(p, p) for p in orients] )
        Q = np.mean(Qs, axis=0) - np.identity(3)/3

        # Nematic order S
        S = np.sqrt( np.tensordot(Q, Q)*1.5)

        for syl in self.sylinders:
            syl.polarOrder = polarOrder
            syl.nematicOrder = S
            syl.nematicTensor = Q

        return S, polarOrder

    def parseSylinderFile(self, sylinderFile):
        self.parseFile(sylinderFile, Sylinder, self.sylinders)

    # member variables are dynamically added by parsing data files
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
            syl = objType()
            syl.end0 = np.array( data.GetPoints().GetPoint(2 * i) )
            syl.end1 = np.array( data.GetPoints().GetPoint(2 * i + 1) )
            xi = syl.end1 - syl.end0 
            syl.orientation = xi/np.sqrt(xi.dot(xi))
            objList.append(syl)

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


if __name__ == '__main__':
    print('not implemented yet')    
