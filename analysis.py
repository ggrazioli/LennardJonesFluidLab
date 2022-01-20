#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli
"""  

import pandas as pd  
import numpy as np
from mathTools import PeriodicVectorCorrection
    
def GetRDF(ensemble, radSteps=20, sampleFrac=.4):
    '''
    Parameters
    ----------
    ensemble : Ensemble
        Ensemble class.
    radSteps : int, optional
        Number of steps in radius value for RDF. The default is 20.
    sampleFrac : float, optional
        Fraction of atoms to select for RDF (<= 1). The default is .2.

    Returns
    -------
    Pandas DataFrame
        Radial distribution function (RDF). 
    '''
    # Returns a radial distribution function for the ensemble ranging 
    # from zero to half the length of the side of the container 
    sampleNum = round(sampleFrac * ensemble.nMol)
    def CountMols(ensemble, step, moID, radius):
        m1 = ensemble.trajectory[step][moID] 
        ct = 0
        for m,mol in enumerate(ensemble.trajectory[step]):
            diff = mol - m1 
            if ensemble.periodic:
                diff = PeriodicVectorCorrection(diff, ensemble.sideLen)
            if ((np.sqrt(sum((diff)**2)) <= radius)  and (moID != m) ):
                ct += 1
        return(ct)
    def GetOneRadius(ensemble, radius, sampleFrac=sampleFrac):
        out = 0. 
        for t in range(0, len(ensemble.trajectory)):
            sampledIndices = np.random.choice(list(range(ensemble.nMol)), sampleNum, replace=False) 
            for m in sampledIndices:
                out += CountMols(ensemble, t, m, radius)
        return (out / (len(ensemble.trajectory) * sampleNum)) 
    RDF = []
    molarDensity = ensemble.nMol / (ensemble.sideLen ** ensemble.dim)
    radList = np.linspace(ensemble.sideLen/9999999, ensemble.sideLen/2., radSteps)
    for idx,r in enumerate(radList):
        if idx % 5 == 0: print("RDF calculation", '{0:.0f}'.format(idx/len(radList)*100), "% complete.")
        nExpect = molarDensity * (ensemble.dim + 1)/3. * np.pi * r ** ensemble.dim
        RDF.append( GetOneRadius(ensemble, r) / nExpect)
    out = zip(np.linspace(0, ensemble.sideLen/2, radSteps), RDF)
    return (pd.DataFrame(out, columns=['r','RDF']))

    
