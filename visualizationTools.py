#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli

"""

import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
    
def makeTrajMovie2D(traj, sideLen, filename = 'LJtraj.gif'):
    fig = plt.figure(figsize=(5,5))
    plt.xlim(0 - .5, sideLen+.5)
    plt.ylim(0 - .5, sideLen+.5)
    graph, = plt.plot([], [], 'o')
    def animate(i):
        x = traj[i][:,0]
        y = traj[i][:,1]
        graph.set_data(x, y)
        return(graph,)
    ani = FuncAnimation(fig, animate, frames=len(traj), blit=True)
    writergif = PillowWriter(fps=30) 
    ani.save(filename, writer=writergif)
    
def plotKEtotals(filepath):
    KEdf = pd.read_csv(filepath, header=0)
    KEdf.plot(x='t', y='KE')
    
def plotPEtotals(filepath):
    PEdf = pd.read_csv(filepath, header=0)
    PEdf.plot(x='t', y='PE')
    
def plotKEandPE(filepathKE, filepathPE):
    KEdf = pd.read_csv(filepathKE, header=0)
    PEdf = pd.read_csv(filepathPE, header=0)
    df = pd.concat([KEdf, PEdf.PE], axis=1)
    df.plot(x='t', y=['KE','PE'])
    
def plotTotalEnergy(filepathKE, filepathPE):
    KEdf = pd.read_csv(filepathKE, header=0)
    PEdf = pd.read_csv(filepathPE, header=0)
    df = pd.concat([KEdf.KE, PEdf.PE], axis=1)
    dfForPlot = pd.concat([KEdf.t, df.sum(axis=1)], axis=1)
    dfForPlot.columns = ['t', 'totEnergy']
    dfForPlot.plot(x='t', y='totEnergy')
    
def plotKEandPEandTotal(filepathKE, filepathPE): 
    KEdf = pd.read_csv(filepathKE, header=0)
    PEdf = pd.read_csv(filepathPE, header=0)
    df = pd.concat([KEdf.KE, PEdf.PE], axis=1)
    dfForPlot = pd.concat([KEdf.t, df, df.sum(axis=1)], axis=1)
    dfForPlot.columns = ['t', 'KE', 'PE', 'totEnergy']
    dfForPlot.plot(x='t', y=['KE','PE', 'totEnergy'])
    
def plotLJpotential(ensemble, ptCount=200):
    x = np.linspace(ensemble.sigma, ensemble.cutoff, ptCount)
    y = []
    for xVal in x:
        sigOverR6 = np.power(ensemble.sigma/xVal, 6)
        sigOverR12 = sigOverR6*sigOverR6        
        y.append( 4 * ensemble.epsilon * (sigOverR12 - sigOverR6))
    df = pd.DataFrame()
    df['distance'] = x
    df['PE'] = y
    df.plot(x = 'distance', y = 'PE')
    
def writeXYZ(ensemble, filename, atomName='Ar'):
    f = open(filename, 'w')
    molNum = str(ensemble.nMol) 
    is2D = (len(ensemble.trajectory[0][0]) == 2)
    for step in ensemble.trajectory:
        f.write(molNum+'\n\n')
        for mol in step:
            f.write(atomName+'\t')
            for dim in mol:
                f.write('%.4f    ' % dim)
            if (is2D): 
                f.write("0.0000")
            f.write('\n')
    f.close()     
    
    