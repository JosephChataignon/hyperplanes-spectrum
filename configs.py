# -*- coding: utf-8 -*-
import numpy as np
import scipy

def getNumberOfRegions( nHyperplanes ):
    return 1 + nHyperplanes + scipy.special.binom(nHyperplanes,2)

def divideRegion( n, region, sideA, sideB, newIndexA, newIndexB ):
    '''
        Divides the region region in 2 other regions, with the new hyperplane n
        passing by sides of index sideA and sideB. newIndexA and newIndexB are 
        the indices of the 2 new regions formed.
    '''
    region1 = [[],[]]
    region2 = [[],[]]
    firstRegion = True
    for k in xrange(len(region[0])):
        if region[0,k] == sideA or region[0,k] == sideB:
            region1[0].append(region[0,k])
            region1[1].append(region[1,k])
            region2[0].append(region[0,k])
            region2[1].append(region[1,k])
            if firstRegion:
                region1[0].append(newIndexB)
                region1[1].append(n)
            else:
                region2[0].append(newIndexA)
                region2[1].append(n)
            firstRegion = not firstRegion # toggle firstRegion value
        else:
            if firstRegion:
                region1[0].append(region[0,k])
                region1[1].append(region[1,k])
            else:
                region2[0].append(region[0,k])
                region2[1].append(region[1,k])
    return region1,region2

def recursiveConfigs( n, config, indexRegion, previousRegion, hyperplanesCrossed=[] ):
    #fin de parcours quand on a traverse tous les hyperplans et que la derniere region est infinie
    newConfigs = []
    region = config[indexRegion]
    
    for i in xrange(len(region[0])):
        if region[1,i] not in hyperplanesCrossed:
            #for each region i neighboring region
            
            if len(hyperplanesCrossed) == n:
                if region[0,i] == -1:
                    #fin
                else:
                    print "!!! !!! Error: blocked by hyperplanes already crossed"
                    return []
            else:
                #nouvelle region traversable: couper la region en 2
                newConfig = config
                region1,region2 = divideRegion( n, region, region[0].index(previousRegion), i, indexRegion, len(config) )
                newConfig[indexRegion] = region1
                newConfig.append(region2)
                #remplacer les references a cette region dans tout le tableau newConfig
                for j in xrange(len(region2[0])):
                    if region2[0,j] != -1:
                        
                #TODO suggestion: utiliser l'ordre (connu) des cotes 
                #et la connaissance du dernier cote traverse pour comment labelliser la region suivante
                
                
                
                hyperplanesCrossed.append(region[1,i])
                #appeler
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
            for k in config:
                departureRegion = config[k]
                if -1 in departureRegion:
                    newConfigs = recursiveConfigs(n,config,indexRegion=k,previousRegion=-1)


























