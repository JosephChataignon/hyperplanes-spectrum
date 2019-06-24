# -*- coding: utf-8 -*-
import numpy as np
import scipy
import copy

def getNumberOfRegions( nHyperplanes ):
    return 1 + nHyperplanes + scipy.special.binom(nHyperplanes,2)

def divideRegion( n, region, sidePreviousRegion, sideNextRegion, previousSecondRegion, newIndexA, newIndexB ):
    '''
        Divides the region region in 2 other regions, with the new hyperplane n
        passing by sides of local index sidePreviousRegion and sideNextRegion.
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
        if (not started) and (k == sidePreviousRegion):
            started = True
            firstRegion = True
            region1[0].append(region[0][k])
            region1[1].append(region[1][k])
        elif started and (k == sideNextRegion):
            firstRegion = False
            region1[0].append(region[0][k])
            region1[1].append(region[1][k])
            region1[0].append(newIndexB)
            region1[1].append(n)
            region2[0].append(newIndexA)
            region2[1].append(n)
            region2[0].append(region[0][k])
            region2[1].append(region[1][k])
        elif started and (k != sideNextRegion) and (k != sideNextRegion):
            if firstRegion:
                region1[0].append(region[0][k])
                region1[1].append(region[1][k])
            else:
                region2[0].append(region[0][k])
                region2[1].append(region[1][k])
        elif started and (k == sidePreviousRegion):
            region2[0].append(previousSecondRegion)
            region2[1].append(region[1][k])
            break
        if started:
            c = c + 1
    print 'region1'
    return region1,region2

def recursiveConfigurations( n, config, indexRegion, previousRegion=-1, previousSecondRegion=-1, hyperplanesCrossed=[] ):
    '''
        Recursive function to generate configurations. n is the index of the 
        hyperplane being added.
    '''
    
    print '\ncall to recursiveConfigs\n n =',n,'\n config =',config,'\n hyperplanesCrossed =',hyperplanesCrossed
    
    # end when the n-th hyperplanes already crossed the n-1 other hyperplanes
    if len(hyperplanesCrossed) >= n :
        print 'meet stop condition, n=',n
        if indexRegion == -1:
            return [config]
        else:
            print "error: n hyperplanes crossed but not arrived at an infinite region"
    else:
        print 'recursion n=',n
        newConfigs = []
        region = config[indexRegion]
        
        for localIndex in xrange(len(region[0])):
            if ( region[1][localIndex] not in hyperplanesCrossed ) and ( region[0][localIndex] != -1 ):
                # for each neighboring region at localIndex reachable (without 
                # re-crossing a hyperplane) and not equal to infinity
                # cut the current region to reach that next one
                
                newConfig = config
                region1,region2 = divideRegion( n, copy.deepcopy(region),
                    region[0].index(previousRegion), localIndex, 
                    previousSecondRegion, indexRegion, len(config) )
                newConfig[indexRegion] = region1
                newConfig.append(region2)
                
                # replace all references to that region in the array newConfig
                for j in xrange(len(region2[0])):
                    # for all surrounding regions of region2 except infinity, region1 and the next region
                    if region2[0][j] not in [ -1 , indexRegion , region[0][localIndex] ]:
                        modifiedRegion = newConfig[region2[0][j]]
                        z = modifiedRegion[0].index(indexRegion)
                        modifiedRegion[0][z] = len(newConfig)-1 # index of region2
                        newConfig[region2[0][j]] = modifiedRegion
                
                # recursively call the function and append its result to newConfigs
                nextHyperplanesCrossed = hyperplanesCrossed + [region[1][localIndex]]
                recursiveConfigs = recursiveConfigurations( n+1, copy.deepcopy(newConfig),
                    region[0][localIndex], indexRegion, 
                    len(newConfig), copy.deepcopy(nextHyperplanesCrossed) )
                print 'recursiveConfigs: ', recursiveConfigs
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
    config  = [region0,region1]
    nextConfigs = [config]
    # Start with a config of 1 hyperplane
    for n in xrange(1,nHyperplanes+1):
        print 'boucle1,n=',n
        configs = nextConfigs
        nextConfigs = []
        for config in configs:
            print 'boucle2,config=',config
            for k in xrange(len(config)):
                print 'boucle3,k=',k
                departureRegion = config[k]
                if -1 in departureRegion[0]:
                    newConfigs = recursiveConfigurations(n,copy.deepcopy(config),indexRegion=k)
                    nextConfigs = nextConfigs + newConfigs
    # TODO :
    # Check in configs if 2 configs differ only by hyperplanes having inverted 
    # departure and arrival regions, or if configs are identical
    # This can be done between adding new hyperplanes so less equivalent configs are generated
    print '\nFinal configs:\n',nextConfigs
    
generateConfigurations(4)

























