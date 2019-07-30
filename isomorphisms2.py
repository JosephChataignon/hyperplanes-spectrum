# -*- coding: utf-8 -*-
# Control program to make sure isomorphisms is working correctly, and that 
# a graph isomorphism is equivalent to identical configurations. 
# In addition to the normal Weinberg's algorithm, we add to the vector a second 
# half detailing which hyperplanes are crossed and in which order.

def getNodesValence(configs):
    """
    Outputs a 2-dimensions array containing, for each configuration, a list of
    the numbers of nodes having the valence (ie, number of branches) 
    corresponding to the index in the list.
    For example, if there are 3 nodes with valence equal to 5 in the graph of 
    the i-th config, then the output at coordinates (i,5) will be 3.
    Applies to a config processed by addInfiniteRegion()
    """
    nodesValence = []
    for config in configs:
        val = [0]*(len(config)+1)
        for region in config:
            val[ len(region[0]) ] += 1
        nodesValence.append(val)
    return nodesValence

def choseNextBranch(branchesVisited,node,previousBranch,reverseOrder):
    """
    Returns the next branch that has not already been used in this direction, 
    at node node and coming from branch previousBranch.
    Default order is clockwise, but if reverseOrder is true the order is
    counter-clockwise.
    """
    started = False
    for k in range( 2*len(branchesVisited[node]) + 1 ):
        if reverseOrder:
            k = -k
        b = k % len(branchesVisited[node])
        if started and (not branchesVisited[node][b]):
            return b
        if b == previousBranch:
            started = True    

def generateCodeVector(config,r,b,reverseOrder=False):
    """
    Given a configuration, a node of index r and a branch of index b on this 
    node, this function follows Weinberg's algorithm to generate the eulerian 
    path and return the code vector of the configuration config for starting 
    branch b on node r.
    Variable reverseOrder indicates whether the order is default or not 
    (clockwise or counterclockwise).
    Applies to a config processed by addInfiniteRegion()
    """
    nodesVisited = [False]*len(config); nodesVisited[r] = True
    branchesVisited = [ [False]*len(region[0]) for region in config]
    path = [r]; pathVertices = []
    initialNode = r; branch = b; terminalNode = config[r][0][b];
    
    while(True):
        reverseBranch = config[terminalNode][0].index(initialNode)
        
        branchesVisited[initialNode][branch] = True
        path.append(terminalNode)
        pathVertices.append(config[initialNode][1][branch])
        
        if nodesVisited[terminalNode] == False: #new node: go to next branch
            nodesVisited[terminalNode] = True
            
            branch = choseNextBranch(branchesVisited,terminalNode,reverseBranch,reverseOrder)
            if branch == None:
                print("Error: no available branch on a new node")
                break
            initialNode = terminalNode
            terminalNode = config[terminalNode][0][branch]
        
        else: #old node
            
            if branchesVisited[terminalNode][reverseBranch] == False: #new branch: make a U-turn
                branchesVisited[terminalNode][reverseBranch] = True
                path.append(initialNode)
                pathVertices.append(config[terminalNode][1][reverseBranch])
                branch = choseNextBranch(branchesVisited,initialNode,branch,reverseOrder)
                if branch == None:
                    print("Error: no available branch while coming from a new branch\n")
                    print("path:\n",path,"\n")
                    break
                terminalNode = config[initialNode][0][branch]
                
            else: #old branch: go to next available branch
                branch = choseNextBranch(branchesVisited,terminalNode,reverseBranch,reverseOrder)
                if branch == None:
                    break
                initialNode = terminalNode
                terminalNode = config[terminalNode][0][branch]
    else:
        print("Error: no break instruction reached in the loop")
    correspondance = []
    vector = []
    vectorVertices = []
    for n in path:
        if not( n in correspondance ):
            correspondance.append(n)
            vector.append( len(correspondance)-1 )
        else:
            vector.append( correspondance.index(n) )
    correspondance = []
    for n in path:
        if not( n in correspondance ):
            correspondance.append(n)
            vectorVertices.append( len(correspondance)-1 )
        else:
            vectorVertices.append( correspondance.index(n) )
    return vector+vectorVertices

def findGraphEdge(config,infiniteIndex):
    """
    Returns the first node of config that is on the border of the graph (ie, 
    that is linked by an edge to the infinite region). infiniteIndex is the 
    index of the infinite region in config.
    Applies to a config processed by addInfiniteRegion()
    """
    for r in range(len(config)):
        for b in config[r][0]:
            if b == infiniteIndex:
                return r

def addInfiniteRegion(configs):
    """
    Adds the infinite region as a node of the graph to each configuration in
    the array configs
    """
    infiniteIndex = len(configs[0])
    # replace reference to the infinite region and delete hyperplanes
    configs = [[[[infiniteIndex if b == -1 else b for b in ri] for ri in r] for r in config] for config in configs]
    # add the new region
    for k in range(len(configs)):
        newRegion = [[],[]]
        r = findGraphEdge(configs[k],infiniteIndex)
        while True:
            localR = (configs[k][r][0].index(infiniteIndex)-1)%len(configs[k][r][0])
            v = configs[k][r][1][localR]
            r = configs[k][r][0][localR]
            if r in newRegion[0]:
                break
            newRegion[0].append(r)
            newRegion[1].append(v)
        configs[k].append(newRegion)
    return configs

def graphBranchFromLabel(config,n):
    """
    Returns indices of the node and branch of the branch corresponding to label
    n in the graph of config (with oriented branches, hence the variable 
    reverseOrder)
    returns -1 as r if n > number of branches
    Applies to a config processed by addInfiniteRegion()
    """
    r = 0; b = -1; reverseOrder = False
    for c in range(n+1):
        b += 1
        if len(config[r][0]) <= b:
            b = 0; r += 1
        if len(config) <= r:
            if reverseOrder == False:
                b = 0; r = 0; reverseOrder = True
            else:
                return -1,-1,False
    return r,b,reverseOrder

def checkIdenticalExistingVectors( newConfigs, vectors, v ):
    """
    Returns True if there exists an already generated vector for one of the 
    newConfigs vectors, that is equal to v.
    Returns None otherwise.
    """
    for i in newConfigs:
        for j in vectors[i]:
            if j == v:
                return True

def checkIdenticalNewVectors( configs, newConfigs, vectors, v ):
    """
    Generates new vectors and returns True if one of them is equal to v.
    Returns None otherwise.
    """
    # to do: add a valence check before generating vectors
    for i in newConfigs:
        while True:
            
            # Get the next vector for config i
            r,b,reverseOrder = graphBranchFromLabel( configs[i],len(vectors[i]) )
            if r == -1:
                break # if all vectors have already been generated for i, get to next new config
            newVect = generateCodeVector(configs[i],r,b,reverseOrder)
            
            # Add it to vectors
            vectors[i].append(newVect)
            # Check if it is equal to v
            if newVect == v:
                return True

def checkValence( k,newConfigs,nodesValence ):
    """
    Returns True if config k's valence is different from any other new config's
    valence, False otherwise.
    """
    for i in newConfigs:
        # check if valences are equal
        if nodesValence[k] == nodesValence[i]:
            return False
    return True

def Weinberg( configs ):
    """
    Applies the Weinberg algorithm to reduce the size of configs by eliminating
    a config when 2 of them are isomorphic.
    """
    # add region -1
    configsInf = addInfiniteRegion(configs)
    
    # generate the vector of nodes valence and the vector of meshes shapes
    nodesValence = getNodesValence(configs)
    
    #meshesShapes = getMeshesShapes(configs)
    # To do: we still need to find an efficient way to get meshes shapes
    
    
    newConfigs = [0] # contains INDICES from configs
    vectors = [ [] for k in range(len(configsInf)) ] # list of the lists of vectors of each graph/config
    
    # for each config k, check if it is isomorphic to a previous config i
    for k in range( 1 , len(configsInf) ):
        
        if checkValence(k,newConfigs,nodesValence):
            newConfigs.append(k)
            
        # if there is at least one graph from newConfigs with same valence as k
        else: 
            v = generateCodeVector(configsInf[k],0,0)
            # check existing vectors compared to v
            if checkIdenticalExistingVectors( newConfigs, vectors, v ):
                continue
                    
            # if no already generated vectors match with v, generate vectors to previously registered configs
            elif checkIdenticalNewVectors(configsInf, newConfigs, vectors, v):
                continue
            else:
                newConfigs.append(k)
                vectors[k] = [v]
    
    # Conversion from indices to actual configurations
    configsReturned = []
    for i in newConfigs:
        configsReturned.append(configs[i])
    
    return configsReturned









## *** Tests area ***

#config 2 regions:
#config = [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]]
#config2 = [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]]

#configs=[config,config2]

#configs = [[[[-1, 1, 11], [-1, 0, 4]], [[0, -1, 3, 12], [0, -1, 1, 4]], [[4, 15, -1], [1, 0, -1]], [[1, -1, 6, 13], [1, -1, 2, 4]], [[7, 5, 2, -1], [2, 0, 1, -1]], [[8, 15, 4], [2, 1, 0]], [[3, -1, 10, 14], [2, -1, 3, 4]], [[11, 8, 4, -1], [3, 0, 2, -1]], [[12, 9, 5, 7], [3, 1, 2, 0]], [[13, 15, 8], [3, 2, 1]], [[6, -1, 15], [3, -1, 4]], [[0, 12, 7, -1], [4, 0, 3, -1]], [[1, 13, 8, 11], [4, 1, 3, 0]], [[3, 14, 9, 12], [4, 2, 3, 1]], [[6, 15, 13], [4, 3, 2]], [[10, -1, 2, 5, 9, 14], [4, -1, 0, 1, 2, 3]]], 
#           [[[-1, 1, 11], [-1, 0, 4]], [[0, -1, 3, 12], [0, -1, 1, 4]], [[4, 15, -1], [1, 0, -1]], [[1, -1, 6, 9, 13], [1, -1, 2, 3, 4]], [[7, 5, 2, -1], [2, 0, 1, -1]], [[8, 15, 4], [2, 1, 0]], [[3, -1, 10], [2, -1, 3]], [[11, 8, 4, -1], [3, 0, 2, -1]], [[12, 14, 5, 7], [3, 1, 2, 0]], [[3, 10, 14], [3, 2, 4]], [[9, 6, -1, 15], [2, 3, -1, 4]], [[0, 12, 7, -1], [4, 0, 3, -1]], [[1, 13, 8, 11], [4, 1, 3, 0]], [[3, 14, 12], [4, 3, 1]], [[9, 15, 8, 13], [4, 2, 1, 3]], [[10, -1, 2, 5, 14], [4, -1, 0, 1, 2]]], 
#           [[[-1, 1, 11], [-1, 0, 4]], [[0, -1, 3, 8, 12], [0, -1, 1, 3, 4]], [[4, 15, -1], [1, 0, -1]], [[1, -1, 6, 9], [1, -1, 2, 3]], [[7, 5, 2, -1], [2, 0, 1, -1]], [[13, 15, 4], [2, 1, 0]], [[3, -1, 10], [2, -1, 3]], [[11, 13, 4, -1], [3, 0, 2, -1]], [[1, 9, 13], [3, 1, 4]], [[8, 3, 10, 14], [1, 3, 2, 4]], [[9, 6, -1, 15], [2, 3, -1, 4]], [[0, 12, 7, -1], [4, 0, 3, -1]], [[1, 13, 11], [4, 3, 0]], [[8, 14, 5, 7, 12], [4, 1, 2, 0, 3]], [[9, 15, 13], [4, 2, 1]], [[10, -1, 2, 5, 14], [4, -1, 0, 1, 2]]],
#           [[[-1, 12, 14], [-1, 0, 3]], [[3, 8, 12], [1, 3, 4]], [[4, 10, -1], [1, 0, -1]], [[-1, 6, 9, 1, 11], [-1, 2, 3, 1, 4]], [[7, 5, 2, -1, 15], [2, 0, 1, -1, 4]], [[8, 10, 4], [2, 1, 0]], [[3, -1, 10], [2, -1, 3]], [[8, 4, 14], [0, 2, 4]], [[1, 9, 5, 7, 13], [3, 1, 2, 0, 4]], [[3, 10, 8], [3, 2, 1]], [[6, -1, 2, 5, 9], [3, -1, 0, 1, 2]], [[3, 12, -1], [4, 1, -1]], [[1, 13, 0, -1, 11], [4, 3, 0, -1, 1]], [[8, 14, 12], [4, 0, 3]], [[7, 15, -1, 0, 13], [4, 2, -1, 3, 0]], [[4, -1, 14], [4, -1, 2]]]]

#print(len(configs))
#z=Weinberg(configs)
#
#print(len(z))












