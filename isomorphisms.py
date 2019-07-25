# -*- coding: utf-8 -*-




# The following 3 functions are from or adapted from :
# github.com/ddatsko/isomorphic_graphs
# however this part of the program is highly inefficient in this case
def next_permutation(perm):
    """
    Precondition: all the elements must be different
    Generates the next permutation by the given one
    Returns None if the permutation is the last
    """
    n = len(perm)
    j = n - 2
    while perm[j] > perm[j + 1]:
        j -= 1
    k = n - 1
    while perm[j] > perm[k]:
        k -= 1
    perm[j], perm[k] = perm[k], perm[j]
    r = n - 1
    s = j + 1
    while r > s:
        perm[r], perm[s] = perm[s], perm[r]
        r -= 1
        s += 1
    if sorted(perm) == perm:
        return None
    return perm

def formatConfig(config):
    """
    config -> dict
    Converts from format [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]] to dictionary
    """
    dct_vertices = {'-1': []}
    for vertex in range(len(config) ):
        dct_vertices[str(vertex)] = []
        for linkedVertex in config[vertex][0]:
            dct_vertices[str(vertex)].append( str(linkedVertex) )
            if linkedVertex == -1:
                dct_vertices['-1'].append( str(vertex) )
    return dct_vertices

def check_isomorphism(graph1, graph2):
    """
    (dict, dict) -> bool
    Check if 2 graphs are isomorphic
    """
    # Check if they are not isomorphic by simple signs(number of edges,
    # number of vertices, degrees of all the vertices
    if len(graph1) != len(graph2):
        return False
    degrees = [0] * len(graph1) * 2
    for i in graph1:
        degrees[len(graph1[i])] += 1
    for i in graph2:
        degrees[len(graph2[i])] -= 1
    for vertice in degrees:
        if vertice != 0:
            return False

    # check by considering all the possible permutations of vertices
    perm = sorted(list(graph2.keys()))
    while perm:
        bad = -1
        change = dict(zip(graph1.keys(), perm))
        for key in graph1:
            if len(graph1[key]) != len(graph2[change[key]]):
                bad = change[key]
                break
            for vertex in graph1[key]:
                if change[vertex] not in graph2[change[key]]:
                    break
            else:
                continue
            break
        else:
            return change
        if bad != -1 and perm.index(bad) != len(perm)-1:
            perm[perm.index(
                bad)+1:] = sorted(perm[perm.index(bad)+1:], reverse=True)
        perm = next_permutation(perm)
    return False

def bruteForce( configs ):
    """
    Function that compares configurations to eliminate the doubles (which graph
    is isomorphic to another configuration graph)
    """
    graphs = [formatConfig(config) for config in configs]
    newConfigs = []
    for i in range(len(graphs) ):
        for j in range(i):
            if check_isomorphism(graphs[i],graphs[j]):
                break
        else:
            newConfigs.append(configs[i])
    return newConfigs



###############################################################################

# This second method is based on the algorithm presented by Louis Weinberg in 
# his paper "A Simple and Efficient Algorithm for Determining Isomorphism of 
# Planar Triply Connected Graphs" from 1965, available at 
# ieeexplore.ieee.org/document/1082573


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
            val[ len(region) ] += 1
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
    branchesVisited = [ [False]*len(region) for region in config]
    path = [r]
    initialNode = r; branch = b; terminalNode = config[r][b];
    # Inside the loop, we mark the current branch (between initialNode and 
    # terminalNode) as visited, as well as terminalNode if it isn't already,
    # we add terminalNode to path, and update the value of initialNode,
    # terminalNode and branch.
    while(True):#Could be changed, length of eulerian path
        reverseBranch = config[terminalNode].index(initialNode)
        
        branchesVisited[initialNode][branch] = True
        path.append(terminalNode)
        
        if nodesVisited[terminalNode] == False: #new node: go to next branch
            nodesVisited[terminalNode] = True
            
            branch = choseNextBranch(branchesVisited,terminalNode,reverseBranch,reverseOrder)
            if branch == None:
                print("Error: no available branch on a new node")
                break
            initialNode = terminalNode
            terminalNode = config[terminalNode][branch]
        
        else: #old node
            
            if branchesVisited[terminalNode][reverseBranch] == False: #new branch: make a U-turn
                branchesVisited[terminalNode][reverseBranch] = True
                path.append(initialNode)
                branch = choseNextBranch(branchesVisited,initialNode,branch,reverseOrder)
                if branch == None:
                    print("Error: no available branch while coming from a new branch\n")
                    print("path:\n",path,"\n")
                    break
                terminalNode = config[initialNode][branch]
                
            else: #old branch: go to next available branch
                branch = choseNextBranch(branchesVisited,terminalNode,reverseBranch,reverseOrder)
                if branch == None:
                    break
                initialNode = terminalNode
                terminalNode = config[terminalNode][branch]
    else:
        print("Error: no break instruction reached in the loop")
    correspondance = []
    vector = []
    for n in path:
        if not( n in correspondance ):
            correspondance.append(n)
            vector.append( len(correspondance)-1 )
        else:
            vector.append( correspondance.index(n) )
    return vector

def generateCodeMatrix(config):
    """
    Generates the full code matrix from a given configuration.
    Applies to a config processed by addInfiniteRegion()
    """
    codeMat = []
    for r in range(len(config)): # r like region
        for b in range(len( config[r] )): # b like branch
            codeMat.append( generateCodeVector(config,r,b) )
            codeMat.append( generateCodeVector(config,r,b,True) )
    # To do: sort the vectors
    return codeMat

def findGraphEdge(config,infiniteIndex):
    """
    Returns the first node of config that is on the border of the graph (ie, 
    that is linked by an edge to the infinite region). infiniteIndex is the 
    index of the infinite region in config.
    Applies to a config processed by addInfiniteRegion()
    """
    for r in range(len(config)):
        for b in config[r]:
            if b == infiniteIndex:
                return r

def addInfiniteRegion(configs):
    """
    Adds the infinite region as a node of the graph to each configuration in
    the array configs
    """
    infiniteIndex = len(configs[0])
    # replace reference to the infinite region and delete hyperplanes
    configs = [[[infiniteIndex if b == -1 else b for b in r[0]] for r in config] for config in configs]
    # add the new region
    for k in range(len(configs)):
        newRegion = []
        r = findGraphEdge(configs[k],infiniteIndex)
        while True:
            localR = (configs[k][r].index(infiniteIndex)-1)%len(configs[k][r])
            r = configs[k][r][localR]
            if r in newRegion:
                break
            newRegion.append(r)
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
        if len(config[r]) <= b:
            b = 0; r += 1
        if len(config) <= r:
            if reverseOrder == False:
                b = 0; r = 0; reverseOrder = True
            else:
                return -1,-1,False
    return r,b,reverseOrder

def checkIdenticalExistingVectors( configs, newConfigs, vectors, v, k ):
    for i in newConfigs:
        for j in vectors[i]:
            if j == v:
                return True

def checkIdenticalNewVectors( configs, newConfigs, vectors, v ):
    for i in newConfigs:
        while True:
            # Get the next vector for config i
            r,b,reverseOrder = graphBranchFromLabel( configs[i],len(vectors[i]) )
            if r == -1:
                break
            newVect = generateCodeVector(configs[i],r,b,reverseOrder)
            # Add it to vectors
            vectors[i].append(newVect)
            # Check if it is equal to v
            if newVect == v:
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
    vectors = [ [] for k in range(len(configsInf)) ] # list of the lists of vectors of each graph
    
    # for each config k, check if it is isomorphic to a previous config i
    for k in range( 1 , len(configsInf) ):
        
        for i in newConfigs:
            # check if valences are equal
            if nodesValence[k] != nodesValence[i]:
                newConfigs.append(k)
                break
            
        else: # if valences are identical, check existing vectors compared to v
            v = generateCodeVector(configsInf[k],0,0)
            if checkIdenticalExistingVectors( configsInf, newConfigs, vectors, v, k ):
                break
                    
            # if no already generated vectors match with v, add vectors to previously registered configs
            if checkIdenticalNewVectors(configsInf, newConfigs, vectors, v):
                break
            else:
                newConfigs.append(k)
                vectors[k] = [v]
    
    # Conversion from indices to actual configurations
    configsReturned = []
    for i in newConfigs:
        configsReturned.append(configs[i])
    
    return configsReturned # return configs with only indices from newConfigs





#config 2 regions:
#config = [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]]
#config2 = [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]]

#config 4 regions:
#config = [ [[-1, 1,6,5], [0]*4], [[-1,2, 0], [0]*3], [[1,-1,3,6], [0]*4], [[-1,4,2], [0]*3], [[-1,5,6,3], [0]*4], [[-1,0,4], [0]*3], [[0,2,4], [0]*3] ]
#config2 = [[[3,1,2], [0]*3], [[-1,5,0,4],[0]*4], [[-1,6,0,5],[0]*4], [[-1,4,0,6],[0]*4], [[-1,1,3],[0]*3], [[-1,2,1],[0]*3], [[-1,3,2],[0]*3] ]

#configs=[config,config2]
#configs = [[[[-1, 1, 11], [-1, 0, 4]], [[0, -1, 3, 12], [0, -1, 1, 4]], [[4, 15, -1], [1, 0, -1]], [[1, -1, 6, 13], [1, -1, 2, 4]], [[7, 5, 2, -1], [2, 0, 1, -1]], [[8, 15, 4], [2, 1, 0]], [[3, -1, 10, 14], [2, -1, 3, 4]], [[11, 8, 4, -1], [3, 0, 2, -1]], [[12, 9, 5, 7], [3, 1, 2, 0]], [[13, 15, 8], [3, 2, 1]], [[6, -1, 15], [3, -1, 4]], [[0, 12, 7, -1], [4, 0, 3, -1]], [[1, 13, 8, 11], [4, 1, 3, 0]], [[3, 14, 9, 12], [4, 2, 3, 1]], [[6, 15, 13], [4, 3, 2]], [[10, -1, 2, 5, 9, 14], [4, -1, 0, 1, 2, 3]]], [[[-1, 1, 11], [-1, 0, 4]], [[0, -1, 3, 12], [0, -1, 1, 4]], [[4, 15, -1], [1, 0, -1]], [[1, -1, 6, 9, 13], [1, -1, 2, 3, 4]], [[7, 5, 2, -1], [2, 0, 1, -1]], [[8, 15, 4], [2, 1, 0]], [[3, -1, 10], [2, -1, 3]], [[11, 8, 4, -1], [3, 0, 2, -1]], [[12, 14, 5, 7], [3, 1, 2, 0]], [[3, 10, 14], [3, 2, 4]], [[9, 6, -1, 15], [2, 3, -1, 4]], [[0, 12, 7, -1], [4, 0, 3, -1]], [[1, 13, 8, 11], [4, 1, 3, 0]], [[3, 14, 12], [4, 3, 1]], [[9, 15, 8, 13], [4, 2, 1, 3]], [[10, -1, 2, 5, 14], [4, -1, 0, 1, 2]]], [[[-1, 1, 11], [-1, 0, 4]], [[0, -1, 3, 8, 12], [0, -1, 1, 3, 4]], [[4, 15, -1], [1, 0, -1]], [[1, -1, 6, 9], [1, -1, 2, 3]], [[7, 5, 2, -1], [2, 0, 1, -1]], [[13, 15, 4], [2, 1, 0]], [[3, -1, 10], [2, -1, 3]], [[11, 13, 4, -1], [3, 0, 2, -1]], [[1, 9, 13], [3, 1, 4]], [[8, 3, 10, 14], [1, 3, 2, 4]], [[9, 6, -1, 15], [2, 3, -1, 4]], [[0, 12, 7, -1], [4, 0, 3, -1]], [[1, 13, 11], [4, 3, 0]], [[8, 14, 5, 7, 12], [4, 1, 2, 0, 3]], [[9, 15, 13], [4, 2, 1]], [[10, -1, 2, 5, 14], [4, -1, 0, 1, 2]]],[[[-1, 1, 9], [-1, 1, 9]], [[-1, 2, 11, 10, 0], [-1, 2, 11, 10, 0]], [[-1, 3, 1], [-1, 3, 1]], [[-1, 4, 12, 11, 2], [-1, 4, 12, 11, 2]], [[-1, 5, 3], [-1, 5, 3]], [[-1, 6, 13, 12, 4], [-1, 6, 13, 12, 4]], [[-1, 7, 5], [-1, 7, 5]], [[-1, 8, 14, 13, 6], [-1, 8, 14, 13, 6]], [[-1, 9, 7], [-1, 9, 7]], [[-1, 0, 10, 14, 8], [-1, 0, 10, 14, 8]], [[1, 15, 9], [1, 15, 9]], [[1, 3, 15], [1, 3, 15]], [[3, 5, 15], [3, 5, 15]], [[5, 7, 15], [5, 7, 15]], [[7, 9, 15], [7, 9, 15]], [[10, 11, 12, 13, 14], [10, 11, 12, 13, 14]]]]

#print(len(configs))
#z=Weinberg(configs)
#print(z)
#print(len(z))












