#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli 
test token
""" 
 
from classes import Ensemble
from dynamics import RunDynamics
import visualizationTools as vis
from analysis import GetRDF
import time

startTime = time.time()

T_init = 1.0 # initial temperature 

molecules = Ensemble(nMol = 64, sideLen=5, dim=3, sigma=1.0, epsilon=1.0, buffer=.75, periodic=True, useLJ=True, temperature=T_init)              
vis.plotLJpotential(molecules)
            
trajectory = RunDynamics(molecules, nSteps=1000, saveFreq=3)
#vis.makeTrajMovie2D(molecules.trajectory, molecules.sideLen) # best for 2D liquid simulation
vis.plotKEandPEandTotal("KE.csv", "PE.csv")
vis.writeXYZ(molecules, "trajectory.xyz")

RDF = GetRDF(molecules, radSteps=50) 
RDF.plot(x='r', y='RDF')

executionTime = (time.time() - startTime)
print('Execution time: ' + str(int(executionTime / 60)) + " min and ", str(int(executionTime) % 60) + " seconds.")






