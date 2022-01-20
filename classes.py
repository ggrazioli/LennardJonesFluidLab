#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli
"""

import numpy as np
from mathTools import PeriodicVectorCorrection
from mathTools import GetRandomUnitVectors 
from mathTools import GetPointsOnSquareGrid, GetPointsOnCubicGrid

class Mol():
    def __init__(self, r, v, a, m):
        self.r = r   # position vector 
        self.v = v   # velocity vector
        self.a = a   # acceleration vector
        self.m = m   # mass

class Ensemble():
    def __init__(self, sideLen=10, dim=3, buffer=.75, nMol=64, temperature=1, sigma=1., epsilon=1., mass=1.0, cutFactor=3, periodic=True, useLJ=True):
        self.bounds = [(0,sideLen) for i in range(dim)] # boundaries of the container
        self.temperature = temperature
        self.sideLen = sideLen   # length of side of container 
        self.dim = dim  # number of dimensions
        self.sigma = sigma     # Lennard-Jones parameter
        self.epsilon = epsilon # Lennard-Jones parameter
        self.nMol = nMol # number of molecules
        self.cutoff = self.sigma * cutFactor # interaction distance cutoff
        self.periodic = periodic # use periodic boundary conditions 
        self.useLJ = useLJ # True turns on Lennard-Jones potential
        if dim == 2:
            positions = GetPointsOnSquareGrid(nMol, sideLen, buffer)
        else:
            positions = GetPointsOnCubicGrid(nMol, sideLen, buffer)
        v = np.zeros(dim) 
        a = np.zeros(dim)
        self.moList = [] # list of Mol() instances
        for p in positions:
            self.moList.append(Mol(p, v, a, mass)) 
        totMass = 0
        for mol in self.moList:
            totMass += mol.m
        self.density = totMass/(sideLen**dim)
        self.trajectory = [np.vstack(positions)]
        self.velocityTraj = [np.vstack([mol.v for mol in self.moList])] # velocities at each step
        self.instPressures = [] # instantaneous pressures calculated from internal virial 
        # See equation  2.53 in Allen and Tildesley for virial expression
        self.pressure = 0 # dynamics function averages instPressures to get pressure 
        self.volume = sideLen ** dim # volume of container
        self.InitializeVelocities()
    
    def InitializeVelocities(self): 
        # Velocities set based on kinetic theory of gases:
        # 1/2 m v^2 = 3/2 kB T (3 for 3 dimensions)
        # We are working with k_B set to 1
        randVecs = GetRandomUnitVectors(self.dim, self.nMol)
        vSum = np.zeros(self.dim)  
        velMags = [np.sqrt(self.dim * self.temperature / mol.m) for mol in self.moList]
        for i,mol in enumerate(self.moList):
            mol.v = velMags[i] * randVecs[i] 
            vSum += mol.v 
        vCentroid = vSum/self.nMol
        for mol in self.moList:
            mol.v -= vCentroid # subtract centroid for stationary center of mass 
                
    def UpdatePeriodicPositions(self):
        for mol in self.moList:
            pos = mol.r
            for d in range(self.dim):
                if pos[d] > self.sideLen:
                    pos[d] = pos[d] % self.sideLen 
                elif pos[d] < 0:
                    pos[d] = self.sideLen + pos[d]

    def UpdateNonPeriodicPositionsAndVelocities(self):
        for mol in self.moList:
            pos = mol.r
            vel = mol.v
            for d in range(self.dim):
                if pos[d] > self.sideLen:
                    pos[d] =  self.sideLen - (pos[d] - self.sideLen)
                    vel[d] = -vel[d]
                elif pos[d] < 0:
                    pos[d] = -pos[d]
                    vel[d] = -vel[d]
    
    def AddLJpotForce(self, mol1, mol2): 
        # Add Forces due to the Lennard-Jones potential  
        m1, m2 = self.moList[mol1], self.moList[mol2] 
        diff = m2.r - m1.r 
        if self.periodic:
            diff = PeriodicVectorCorrection(diff, self.sideLen)
        r = np.sqrt(sum(diff**2)) 
        sig6 = self.sigma**6 
        sig12 = sig6*sig6 
        r7 = r**7 
        r13 = r7*r7/r 
        LJforce = 48 * self.epsilon * (sig12/r13 - .5 * sig6/r7) 
        m1.a -= LJforce * (diff/r) / m1.m # r is used to normalize the diff vector
        m2.a += LJforce * (diff/r) / m2.m # r is used to normalize the diff vector
              
    def UpdateAccelerationLJ(self):
        # Accelerations due to Lennard-Jones interactions 
        rCutSq = self.cutoff * self.cutoff 
        virSum = 0 # sum for virial treatment of pressure
        vv = 0 # velocity dot products for kinetic treatment of temperature         
        masses = [self.moList[i].m for i in range(0, len(self.moList))]
        for mol in self.moList:
            mol.a = np.zeros(self.dim)
        for i in range(self.nMol):
            vv += np.dot(self.moList[i].v, self.moList[i].v) * masses[i]
            for j in range(self.nMol):
                m1, m2 = self.moList[i], self.moList[j]
                diff = m2.r - m1.r
                if self.periodic:
                    diff = PeriodicVectorCorrection(diff, self.sideLen)
                distSq = sum(diff**2)
                if( (distSq < rCutSq) and (i < j)):
                    self.AddLJpotForce(i, j) 
                    virSum += np.dot(diff, self.moList[i].a)
        self.temperature = vv/(self.nMol * self.dim)
        V, T, D, rho = self.volume, self.temperature, self.dim, self.density 
        P = rho*T + 1/D * virSum/V # See equation 2.55 in Allen and Tildesley
        self.instPressures.append(P)
        
    def UpdatesForNoLJ(self):
        T, rho = self.temperature, self.density 
        P = rho*T # ideal gas for k_B = 1
        self.instPressures.append(P)
        vv = 0 # velocity dot products for kinetic treatment of temperature         
        masses = [self.moList[i].m for i in range(0, len(self.moList))]
        for i in range(self.nMol):
            vv += np.dot(self.moList[i].v, self.moList[i].v) * masses[i]
        self.temperature = vv/(self.nMol * self.dim)
    
    def onePairLJenergy(self, mol1, mol2): 
        # returns Lennard-Jones potential energy between one pair of molecules 
        m1, m2 = self.moList[mol1], self.moList[mol2] 
        diff = m2.r - m1.r 
        if self.periodic:
            diff = PeriodicVectorCorrection(diff, self.sideLen)
        r = np.sqrt(sum(diff**2)) 
        sigOverR6 = np.power(self.sigma/r, 6)
        sigOverR12 = sigOverR6*sigOverR6        
        return (4 * self.epsilon * (sigOverR12 - sigOverR6))
    
    def getPE(self):
        # get total potential energy 
        totPE = 0
        if self.useLJ:
            rCutSq = self.cutoff * self.cutoff  
            for i in range(0, self.nMol):
                for j in range(0, self.nMol):
                    m1, m2 = self.moList[i], self.moList[j]
                    diff = m2.r - m1.r
                    if self.periodic:
                        diff = PeriodicVectorCorrection(diff, self.sideLen)
                    distSq = sum(diff**2)
                    if( (distSq < rCutSq) and (i < j)):
                        totPE += self.onePairLJenergy(i, j) 
        return(totPE) 
                    
    def getKE(self):
        # get total kinetic energy
        vvTot = 0
        for mol in self.moList:
            vvTot += .5 * mol.m * sum(mol.v**2) #1/2 mv^2
        return(vvTot)
    