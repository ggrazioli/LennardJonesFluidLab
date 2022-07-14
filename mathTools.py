#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli
test comment, David, Adam, Heekun
"""

import numpy as np

def GetRandomUnitVectors(dim, count):
    randos = np.random.normal(0, 1, (dim+1)*count).reshape((count, dim+1))
    norms = [np.sqrt(sum(randos[r]**2)) for r in range(len(randos))]
    raw = [randos[i][0:dim]/norms[i] for i in range(len(norms))]
    outNorms = [np.sqrt(sum(raw[r]**2)) for r in range(len(raw))]
    out = [raw[i]/outNorms[i] for i in range(len(outNorms))]
    return(out)
    
def PeriodicVectorCorrection(vec, sideLen): 
    out = [v for v in vec]
    for i in range(len(vec)):
        if vec[i] >= .5*sideLen:
            out[i] -= sideLen
        elif vec[i] < -.5*sideLen:
            out[i] += sideLen
    return (np.array(out))
 
def getLargerPerfectSquare(num):
    i = 0
    diffSq = 9999999
    oldDiffSq = 9999999999
    while diffSq < oldDiffSq:
        oldDiffSq = diffSq
        i += 1
        diffSq = (num - i*i) * (num - i*i)
    if (i-1)*(i-1) >= num:
        out = (i-1)
    else:
        out = i
    return (out)

def getLargerPerfectCube(num):
    i = 0
    diffSq = 9999999
    oldDiffSq = 9999999999
    while diffSq < oldDiffSq:
        oldDiffSq = diffSq
        i += 1
        diffSq = (num - i*i*i) * (num - i*i*i) 
    if (i-1)*(i-1)*(i-1) >= num:
        out = (i-1)
    else:
        out = i
    return (out)

def GetPointsOnSquareGrid(num, sideLen, buffer):
    sideCount = getLargerPerfectSquare(num)
    sideGrid = np.linspace(0+buffer, sideLen-buffer, sideCount)
    fullGrid = []
    for i in sideGrid:
        for j in sideGrid:
            fullGrid.append([i,j])
    return (np.array(fullGrid[0:num]))

def GetPointsOnCubicGrid(num, sideLen, buffer):
    sideCount = getLargerPerfectCube(num)
    sideGrid = np.linspace(0+buffer, sideLen-buffer, sideCount)
    fullGrid = []
    for i in sideGrid:
        for j in sideGrid:
            for k in sideGrid:
                fullGrid.append([i,j,k])
    return (np.array(fullGrid[0:num]))


