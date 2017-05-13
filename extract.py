import numpy as np
import networkx as nx
import utils


def allSimplePathkHop(graph, source, target, k):
    ''' return all the k hops simple paths from nodes in the source list to the target list.
    Parameters
    ----------
    graph : NetworkX graph

    source : node list
       Starting node for path

    target : node list
       Ending node for path. If provided only predecessors between
       source and target are returned

    cutoff : integer, optional
        Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    pred : list of paths
    '''
    hop_k = []
    for node_s in source:
        for node_t in target:
            path_kHop = list(nx.all_simple_paths(graph, node_s, node_t, cutoff=k))
            hop_k = hop_k + path_kHop

    hop_k = [item for item in hop_k if len(item)==k+1] # only keep the paths that their length==k+1
    return hop_k


def get_edges(pathList,k):
    '''extract the all the edges from pathList'''
    edge_layer_all = []
    for i in range(len(pathList)):
        for j in range(k):
            edge_layer_all.append(tuple([np.array(pathList)[i,j],np.array(pathList)[i,j+1]]))
    return edge_layer_all


def extract_node_edges_from_file(fpath1, fpath2, graph,nodes_inv_dict, k):
    '''output the results of subgraph analysis, including the edges, nodes in different layer of the subgraph.
    k: length of path'''
    gene_list1 = utils.readGeneFromFile(fpath1)
    gene_list2 = utils.readGeneFromFile(fpath2)

    # invert the nodes mapping
    # nodes_inv_dict = dict((val, ke) for ke, val in nodes_dict.items())
    key = list(nodes_inv_dict.keys())  # len(key)=1886, load the nodes of immune net

    list1_inter_immune = utils.geneListInter(gene_list1, key)
    list2_inter_immune = utils.geneListInter(gene_list2, key)
    # list1_inter_immune_id = utils.name2id(list1_inter_immune, nodes_inv_dict)
    # list2_inter_immune_id = utils.name2id(list2_inter_immune, nodes_inv_dict)

    pathList = allSimplePathkHop(graph, list1_inter_immune, list2_inter_immune, k)
    edge_layer_all = get_edges(pathList, k)

    if pathList:
        if k==1:
            return list(set(np.array(pathList)[:,0])),list(set(np.array(pathList)[:,1])),edge_layer_all
        elif k==2:
            return list(set(np.array(pathList)[:,0])),list(set(np.array(pathList)[:,1])),\
            list(set(np.array(pathList)[:,2])),edge_layer_all
        elif k==3:
            return list(set(np.array(pathList)[:,0])),list(set(np.array(pathList)[:,1])),\
            list(set(np.array(pathList)[:,2])),list(set(np.array(pathList)[:,3])),edge_layer_all
    else:
        print("There are no interaction between the two groups!!!")
        return