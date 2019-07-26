# -*- coding: utf-8 -*-
# import os
#print os.path.dirname(__file__)
#currentDirectory = os.path.dirname(os.getcwd())
#os.chdir(currentDirectory)
import scipy
import copy

import isomorphisms as isom





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
    while c <= len(region[0]):
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
        elif started and (k != sidePreviousRegion) and (k != sideNextRegion):
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
    return region1,region2

def recursiveConfigurations( n, config, indexRegion, previousRegion=-1, previousSecondRegion=-1, hyperplanesCrossed=[] ):
    '''
        Recursive function to generate configurations. n is the index of the 
        hyperplane being added.
    '''
    
    region = config[indexRegion]
    
    # End when the n-th hyperplanes already crossed the n-1 other hyperplanes
    if len(hyperplanesCrossed) >= n :
        indexNext = region[0].index(-1)
        newConfig = copy.deepcopy(config)
        region1,region2 = divideRegion( n, region=copy.deepcopy(region),
            sidePreviousRegion=region[0].index(previousRegion), sideNextRegion=indexNext, 
            previousSecondRegion=previousSecondRegion, newIndexA=indexRegion, newIndexB=len(config) )
        newConfig[indexRegion] = region1
        newConfig.append(region2)
        # replace all references to that region in the array newConfig
        for j in range(len(region2[0])):
            # for all surrounding regions of region2 except infinity, region1 and the next region
            if region2[0][j] not in [ -1 , indexRegion , region[0][indexNext] ]:
                modifiedRegion = newConfig[region2[0][j]]
                z = modifiedRegion[0].index(indexRegion)
                modifiedRegion[0][z] = len(newConfig)-1 # index of region2
                newConfig[region2[0][j]] = modifiedRegion
        return [newConfig]
        
        
        
    else:
        newConfigs = []
        
        for indexNext in range(len(region[0])):
            if ( region[1][indexNext] not in hyperplanesCrossed ) and ( region[0][indexNext] != -1 ):
                # for each neighboring region at indexNext reachable (without 
                # re-crossing a hyperplane) and not equal to infinity
                # cut the current region to reach that next one
                
                newConfig = copy.deepcopy(config)
                region1,region2 = divideRegion( n, region=copy.deepcopy(region),
                    sidePreviousRegion=region[0].index(previousRegion), sideNextRegion=indexNext, 
                    previousSecondRegion=previousSecondRegion, newIndexA=indexRegion, newIndexB=len(config) )
                newConfig[indexRegion] = region1
                newConfig.append(region2)
                
                # replace all references to that region in the array newConfig
                for j in range(len(region2[0])):
                    # for all surrounding regions of region2 except infinity, region1 and the next region
                    if region2[0][j] not in [ -1 , indexRegion , region[0][indexNext] ]:
                        modifiedRegion = newConfig[region2[0][j]]
                        z = modifiedRegion[0].index(indexRegion)
                        modifiedRegion[0][z] = len(newConfig)-1 # index of region2
                        newConfig[region2[0][j]] = modifiedRegion
                
                # recursively call the function and append its result to newConfigs
                nextHyperplanesCrossed = hyperplanesCrossed + [region[1][indexNext]]
                recursiveConfigs = recursiveConfigurations( n, copy.deepcopy(newConfig), 
                    region[0][indexNext], indexRegion,  
                    len(newConfig)-1, copy.deepcopy(nextHyperplanesCrossed) )
                newConfigs = newConfigs + recursiveConfigs
        return newConfigs




def eliminateDoubles( configs ):
    # return isom.bruteForce(configs)
    return isom.Weinberg(copy.deepcopy(configs))

def generateConfigurations( nHyperplanes ):
    '''
        Generates all the possible configurations of nHyperplanes hyperplanes in
        2 dimensions.
        Configs contains all the configurations. A configuration is modelised by 
        a list of regions. A region is modelised by an ordered list of 
        surrounding regions (represented by their index in the array, -1 
        representing infinity for unbounded regions) and an ordered list of 
        hyperplanes separating these regions (-1 for unbounded regions).
        The lists of surrounding regions and hyperplanes are in clockwise order.
    '''
    # Start with a config of 1 hyperplane
    region0 = [[-1,1],[-1,0]]
    region1 = [[-1,0],[-1,0]]
    config  = [region0,region1]
    nextConfigs = [config]
    # Loop adding hyperplanes
    for n in range(1,nHyperplanes):
        #print('\n\nboucle 1, n=',n)
        configs, nextConfigs = nextConfigs, []
        # Iterate through configurations
        for config in configs:
            #print('\nboucle 2, config=',config)
            # Pass by all "departure points" possible for the new hyperplane
            for departFrom in range(len(config)):
                if -1 in config[departFrom][0]: # if the region is on the edge of infinity
                    #print('boucle 3, departFrom:',departFrom)
                    newConfigs = recursiveConfigurations(n,copy.deepcopy(config),indexRegion=departFrom)
                    nextConfigs = nextConfigs + newConfigs
#        if n==4:
#            print("\nn = ",n,"nextConfigs avant elimination doublons:",len(nextConfigs),"configs")
        
        nextConfigs = eliminateDoubles(nextConfigs)
        print(n+1,'hyperplanes: ',len(nextConfigs),'configuration(s)')
    return nextConfigs
    
    
    

blob=[[[-1, 12, 14], [-1, 0, 3]], [[3, 8, 12], [1, 3, 4]], [[4, 10, -1], [1, 0, -1]], [[-1, 6, 9, 1, 11], [-1, 2, 3, 1, 4]], [[7, 5, 2, -1, 15], [2, 0, 1, -1, 4]], [[8, 10, 4], [2, 1, 0]], [[3, -1, 10], [2, -1, 3]], [[8, 4, 14], [0, 2, 4]], [[1, 9, 5, 7, 13], [3, 1, 2, 0, 4]], [[3, 10, 8], [3, 2, 1]], [[6, -1, 2, 5, 9], [3, -1, 0, 1, 2]], [[3, 12, -1], [4, 1, -1]], [[1, 13, 0, -1, 11], [4, 3, 0, -1, 1]], [[8, 14, 12], [4, 0, 3]], [[7, 15, -1, 0, 13], [4, 2, -1, 3, 0]], [[4, -1, 14], [4, -1, 2]]]

numberOfHyperplanes = 5
z=generateConfigurations(numberOfHyperplanes)
print('final result:',len(z),'configuration(s) for',numberOfHyperplanes,'hyperplanes\n')












