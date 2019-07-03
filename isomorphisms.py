# -*- coding: utf-8 -*-




# The following 3 functions are from or adapted from :
# github.com/ddatsko/isomorphic_graphs
# however this part of the program is highly inefficient in this case
def next_permutation(perm):
    """
    (lst) -> lst
    Precondition: all the elements must be different
    Generates the next permutation by the given one
    Return None if the permutation is the last
    >>> next_permutation([1, 2, 5, 4, 3])
    [1, 3, 2, 4, 5]
    >>> next_permutation([3, 2, 1])
    None
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
    Convert the line of format [[[-1, 1], [-1, 0]], [[-1, 0], [-1, 0]]] to dictionary
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
    graphs = [formatConfig(config) for config in configs]
    newConfigs = []
    for i in range(len(graphs) ):
        for j in range(i):
            if check_isomorphism(graphs[i],graphs[j]):
                break
        else:
            newConfigs.append(configs[i])
    return newConfigs










def Weinberg( configs ):
    # generate the vector of nodes valence and the vector of meshes shapes
    
    # If they have the same characteristics, generate the code matrix
    


















