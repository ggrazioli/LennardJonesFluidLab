#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli
"""  
from numpy import vstack, mean

def takeOneStep(ensemble, dt=.005): #dt was .0005
    '''
    Parameters 
    ----------
    ensemble : Ensemble
        Ensemble class from classes.py.
    dt : float, optional
        time step. The default is .005.

    Returns nothing, but updates the ensemble one time step later.
    ----------
    '''
    oldAccels = [ensemble.moList[i].a.copy() for i in range(0, len(ensemble.moList))]
    for mol in ensemble.moList: #update positions (non-periodic)
        mol.r = mol.r + mol.v*dt + .5*mol.a*dt*dt
    if ensemble.periodic:
        ensemble.UpdatePeriodicPositions()
    else:
        ensemble.UpdateNonPeriodicPositionsAndVelocities()
    if ensemble.useLJ: 
        ensemble.UpdateAccelerationLJ() 
    else:
        ensemble.UpdatesForNoLJ()
    for i,mol in enumerate(ensemble.moList): #update velocities
        mol.v = mol.v + .5*(mol.a + oldAccels[i])*dt 
    
def RunDynamics(ensemble, nSteps = 1000, dt=.005, saveFreq=5): 
    '''
    Parameters
    ----------
    ensemble : Ensemble
        Ensemble class from classes.py..
    nSteps : int, optional
        Number of time steps. The default is 10.
    saveFreq: int, optional
        Interval of time steps at which to save trajectory information
    dt : float, optional
        time step. The default is .005.  

    Returns nothing, but as it runs dynamics, the trajectory is stored in
    ensemble.trajectory, and total kinetic and potential energy is outputted
    to CSV files. 
    -------
    '''
    KE = []
    PE = []
    for t in range(0, nSteps):
        if (t % saveFreq == 0):
            positions = [ensemble.moList[i].r for i in range(0, len(ensemble.moList))]
            velocities = [ensemble.moList[i].v for i in range(0, len(ensemble.moList))]
            ensemble.trajectory.append(vstack(positions))
            ensemble.velocityTraj.append(vstack(velocities))
            KinEn = ensemble.getKE()
            KE.append(KinEn)
            PE.append(ensemble.getPE())
        
        takeOneStep(ensemble, dt=dt)
        ensemble.pressure = mean(ensemble.instPressures) 
        if (t % 100 == 0): print("Simulation step",t,"complete")
    
    with open("KE.csv", "w") as KEfile:
        KEfile.write("t,KE\n") 
        for idx,k in enumerate(KE):
            KEfile.write(str(idx * saveFreq) + ", " + str(k)+"\n")
    with open("PE.csv", "w") as PEfile:
        PEfile.write("t,PE\n")
        for idx,p in enumerate(PE):
            PEfile.write(str(idx * saveFreq) + ", " + str(p)+"\n")
    
