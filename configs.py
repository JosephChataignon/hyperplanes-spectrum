# -*- coding: utf-8 -*-
import numpy as np
import scipy

def getNumberOfRegions( nHyperplanes ):
    return 1 + nHyperplanes + scipy.special.binom(nHyperplanes,2)

def divideRegion( n, region, sidePreviousRegion, sideNextRegion, newIndexA, newIndexB ):
    '''
        Divides the region region in 2 other regions, with the new hyperplane n
        passing by sides of index sidePreviousRegion and sideNextRegion.
        newIndexA and newIndexB are the indices of the 2 new regions formed.
        region1 is the region which perimeter starts after (in clockwise order) 
        sidePreviousRegion and finishing by the hyperplane n.
    '''
    region1 = [[],[]]
    region2 = [[],[]]
    started = False
    c = 0 #counting the number of regions processed
    k = 0
    while c < len(region[0]):            
        k = (k+1) % len(region[0])
        if (not started) and (region[0,k] == sidePreviousRegion):
            started = True
            firstRegion = True
            region1[0].append(region[0,k])
            region1[1].append(region[1,k])
        elif started and (region[0,k] == sideNextRegion):
            firstRegion = False
            
        elif started and (region[0,k] != sideNextRegion) and (region[0,k] != sideNextRegion):
            if firstRegion:
                region1[0].append(region[0,k])
                region1[1].append(region[1,k])
            else:
                region2[0].append(region[0,k])
                region2[1].append(region[1,k])
        elif started and (region[0,k] == sidePreviousRegion):
            
            
        if started:
            c = c + 1
    return region1,region2

def recursiveConfigurations( n, config, indexRegion, previousRegion=-1, previousSecondRegion=-1, hyperplanesCrossed=[] ):
    '''
        Recursive function to generate configurations.
        On the parameters, be careful about indices: some are global indices 
        (of the "config" array), and some are indices of the "region" array
    '''
    # end when the n-th hyperplanes already crossed the n-1 other hyperplanes
    if len(hyperplanesCrossed) >= n-1 and indexRegion == -1:
        return config
    else:
        newConfigs = []
        region = config[indexRegion]
        
        for localIndex in xrange(len(region[0])):
            if region[1,localIndex] not in hyperplanesCrossed:
                # for each neighboring region at localIndex reachable (without re-crossing a hyperplane)
                # cut the current region to reach that next one
                
                newConfig = config
                region1,region2 = divideRegion( n, region, region[0].index(previousRegion), localIndex, indexRegion, len(config) )
                newConfig[indexRegion] = region1
                newConfig.append(region2)
                
                
                #TODO remplacer les references a cette region dans tout le tableau newConfig                        
                #TODO suggestion: utiliser l'ordre (connu) des cotes et la 
                #TODO connaissance du dernier cote traverse pour comment labelliser la region suivante
                for j in xrange(len(region2[0])):
                    if region2[0,j] != -1:
                        
                
                
                
                # recursively call the function and append its result to newConfigs
                nextHyperplanesCrossed = hyperplanesCrossed
                nextHyperplanesCrossed.append(region[1,localIndex])
                recursiveConfigs = recursiveConfigurations( n+1, newConfig, region[0,localIndex],
                    indexRegion, len(newConfig), nextHyperplanesCrossed )
                newConfigs = newConfigs + recursiveConfigs
        return newConfigs

def generateConfigurations( nHyperplanes ):
    '''
        Generates all the possible configurations of nHyperplanes in 2
        dimensions.
        Configs contains all the configurations. A configuration is modelised by 
        a list of regions. A region is modelised by an ordered list of 
        surrounding regions (represented by their index in the array, -1 
        representing infinity for unbounded regions) and an ordered list of 
        hyperplanes separating these regions (-1 for unbounded regions).
        The lists of surrounding regions and hyperplanes are in clockwise order.
    '''
    region0 = [[-1,1],[-1,0]]
    region1 = [[-1,0],[-1,0]]
    config = [region0,region1]
    configs = [config] # Start with a config of 1 hyperplane
    for n in xrange(1,nHyperplanes):
        for config in configs:
            # can check if 2 configs differ only by hyperplanes having inverted departure and arrival regions
            # check for identical configs
            for k in config:
                departureRegion = config[k]
                if -1 in departureRegion:
                    newConfigs = recursiveConfigurations(n,config,indexRegion=k,previousRegion=-1)


























