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
    the i-th config, then the output at coordinates (i,5) will be 3
    """
    nodesValence = []
    for config in configs:
        val = [0]*(len(config)+1)
        for region in config:
            val[ len(region[0]) ] += 1
        nodesValence.append(val)
    return nodesValence

def choseNextBranch(branchesVisited,node,previousNode,reverseOrder):
    """
    Returns the next branch that has not already been used in this direction, 
    at node node and coming from node previousNode.
    Default order is clockwise, but if reverseOrder is true the order is
    counter-clockwise.
    """
    started = False
    for k in range(3*len(branchesVisited[node])):
        if reverseOrder:
            k = -k
        b = k % len(branchesVisited[node])
        if started and (not branchesVisited[node][b]):
            return b
        if b == previousNode:
            started = True    

def generateCodeVector(config,r,b,reverseOrder=False):
    """
    Given a configuration, a node of 
    """
    nodesVisited = [False]*len(config)
    branchesVisited = [ [False]*len(region[0]) for region in config]
    path = [r]
    initialNode = r; branch = b; terminalNode = config[r][0][b];
    # Inside the loop, we mark the current branch (between initialNode and terminalNode) as visited,
    # as well as terminalNode,  add terminalNode to path, and update the value of terminalNode and branch
    for c in range(10000):#To be changed, length of eulerian path
        reverseBranch = config[terminalNode][0].index(initialNode)
        if nodesVisited[terminalNode] == False: #new node
            nodesVisited[terminalNode] = True
            branchesVisited[initialNode][branch] = True
            path.append(terminalNode)
            initialNode = terminalNode
            branch = choseNextBranch(branchesVisited,terminalNode,reverseBranch,reverseOrder)
            if branch == None:
                break
            terminalNode = config[terminalNode][0][branch]
        else: #old node
            if branchesVisited[terminalNode][reverseBranch] == False: #new branch
                branchesVisited[terminalNode][reverseBranch] = True
                branchesVisited[initialNode][branch] = True
                path.append(terminalNode)
                path.append(initialNode)
                initialNode, terminalNode = terminalNode, initialNode
            else: #old branch
                branchesVisited[initialNode][branch] = True
                path.append(terminalNode)
                initialNode = terminalNode
                branch = choseNextBranch(branchesVisited,terminalNode,reverseBranch,reverseOrder)
                if branch == None:
                    break
                terminalNode = config[terminalNode][0][branch]
    else:
        print("No breakpoint reached in the loop")
    correspondance = []
    vector = []
    for n in path:
        if not( n in correspondance ):
            correspondance.append(n)
            vector.append( len(correspondance)-1 )
        else:
            vector.append( correspondance.index(n) )
    #to do: check that correspondance has the same length as config?
    print('path',path)
    return vector

def generateCodeMatrix(config):
    codeMat = []
    for r in range(len(config)): # r like region
        for b in range(len( config[r][0] )): # b like branch
            codeMat.append( generateCodeVector(config,r,b) )
            codeMat.append( generateCodeVector(config,r,b,True) )
    # To do: sort the vectors
    return codeMat

def findGraphEdge(config,infiniteIndex):
    for r in range(len(config)):
        for b in config[r][0]:
            if b == infiniteIndex:
                return r

def addInfiniteRegion(configs):
    infiniteIndex = len(configs[0])
    # replace reference to region - can be used to easily remove hyperplanes
    configs = [[[[infiniteIndex if i == -1 else i for i in j] for j in k] for k in l] for l in configs]
    # add new region
    for k in range(len(configs)):
        newRegion = []
        r = findGraphEdge(configs[k],infiniteIndex)
        while True:
            localR = (configs[k][r][0].index(infiniteIndex)-1)%len(configs[k][r][0])
            r = configs[k][r][0][localR]
            if r in newRegion:
                break
            newRegion.append(r)
        configs[k].append([newRegion,newRegion])
    return configs

def Weinberg( configs ):
    # generate the vector of nodes valence and the vector of meshes shapes
    nodesValence = getNodesValence(configs)
    #meshesShapes = getMeshesShapes(configs) # We still need to find an efficient way to get meshes shapes
    # add region -1
    configs = addInfiniteRegion(configs)
    # for those who have the same characteristics, generate the code matrix
    newConfigs = 











config = [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]]
configs = addInfiniteRegion([config])
print(configs[0])
print( 'vector code', generateCodeVector(configs[0],0,0) )














