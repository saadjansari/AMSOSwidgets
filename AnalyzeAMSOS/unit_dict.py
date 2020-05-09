#!/usr/bin/env python3

import sys, os, pdb

'''
Name: unit_dict.py
Description: Conversion dictionary for AMSOS parameters 
Input: None
Output: None
'''

class UnitDict():
    def __init__(self):
        self.ud ={
            "time" : '(sec)', 
            "timeSnap" : '(sec)', 
            "dt" : '(sec)', 
            "viscosity": '(pN.s.${\mum}^{-2}$)',
            "KBT": '(pN.$\mu$m)',
            "distance" : '($\mu$m)', 
            "nematicOrder" : '', 
            "polarOrder" : '', 
            "end0" : '($\mu$m)',
            "end1" : '($\mu$m)',
            "endLabel0" : '',
            "endLabel1" : '',
            "force" : '(pN)',
            "forceBilateral" : '(pN)',
            "forceCollision" : '(pN)',
            "forceNonBrown" : '(pN)',
            "gid" : '',
            "group" : '',
            "length" : '($\mu$m)',
            "lengthCollision" : '($\mu$m)',
            "omega" : '',
            "omegaBilateral" : '',
            "omegaBrown" : '',
            "omegaCollision" : '',
            "omegaNonBrown" : '',
            "radius" : '($\mu$m)',
            "radiusCollision" : '($\mu$m)',
            "torque" : '(pN.$\mu$m)',
            "torqueBilateral" : '(pN.$\mu$m)',
            "torqueCollision" : '(pN.$\mu$m)',
            "torqueNonBrown" : '(pN.$\mu$m)',
            "vel" : '($\mu$m/sec)',
            "velBilateral" : '($\mu$m/sec)',
            "velBrown" : '($\mu$m/sec)',
            "velCollision" : '($\mu$m/sec)',
            "velNonBrown" : '($\mu$m/sec)',
            "xnorm" : '($\mu$m)',
            "znorm" : '($\mu$m)',
        }
        # 'end0', 'end1', 'endLabel0', 'endLabel1', 'force', 'forceBilateral', 'forceCollision', 'forceNonBrown', 'gid', 'group', 'length', 'lengthCollision', 'omega', 'omegaBilateral', 'omegaBrown', 'omegaCollision', 'omegaNonBrown', 'radius', 'radiusCollision', 'torque', 'torqueBilateral', 'torqueCollision', 'torqueNonBrown', 'vel', 'velBilateral', 'velBrown', 'velCollision', 'velNonBrown', 'xnorm', 'znorm'
    def __getitem__(self, k):
        if k in self.ud:
            return self.ud[k]
        else:
            return ('',1.0, float)

##########################################
if __name__ == "__main__":
    print("Not implemented yet")

