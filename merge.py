import utils


def node_inter_g(dic, gNodeList):
    '''get the intersection between pathway and node_layer. return the intersection as dict, 
    key:pathway name, value:common genelist.
    '''
    mergeNode = []
    for k, v in dic.items():
        intersec = utils.geneListInter(gNodeList, list(v))
        if len(intersec) > 3:
            mergeNode.append((k, intersec))
    return dict(mergeNode)


def node_both_inter_g(dic, list1, list2):
    '''return the nodes from both dic&list1, and dic&list2'''
    mergeNode = []
    for k, v in dic.items():
        intersec1 = utils.geneListInter(list1, list(v))
        intersec2 = utils.geneListInter(list2, list(v))
        if intersec1 and intersec2:
            intersec = intersec1 + intersec2
            if len(intersec) > 3:
                mergeNode.append((k, intersec))
    return dict(mergeNode)


def noninter_patition(sorted_funInterNode):
    '''find the non-intersected partitions of nodes set to merge'''
    geneNode = set(sorted_funInterNode[0][1])
    noninterList = [sorted_funInterNode[0]]
    for item in sorted_funInterNode:
        if not (set(item[1]) & set(geneNode)):
            geneNode = geneNode.union(item[1])
            noninterList.append(item)
    return noninterList


def merge_nodes(G, nodes, new_node, attr_dict=None, **attr):
    """
    Merges the selected `nodes` of the graph G into one `new_node`,
    meaning that all the edges that pointed to or from one of these
    `nodes` will point to or from the `new_node`.
    attr_dict and **attr are defined as in `G.add_node`.
    """

    G.add_node(new_node, attr_dict, **attr)  # Add the 'merged' node

    for n1, n2, data in G.edges(data=True):
        # For all edges related to one of the nodes to merge,
        # make an edge going to or coming from the `new gene`.
        if n1 in nodes:
            G.add_edge(new_node, n2, data)
        elif n2 in nodes:
            G.add_edge(n1, new_node, data)


# for n in nodes: # remove the merged nodes
#         G.remove_node(n)

def copy_merge_nodes(G, G_copy, G_cle, G_ple, nodes, new_node, attr_dict=None, **attr):
    """
    Merge the nodes to new_node, but copy this merge to a new graph G_copy
    attr_dict and **attr are defined as in `G.add_node`.
    """

    #     G_copy.add_node(new_node, attr_dict, **attr) # Add the 'merged' node

    for n1, n2, data in G.edges(data=True):
        # For all edges related to one of the nodes to merge,
        # make an edge going to or coming from the `new gene`.
        if n1 in nodes and n2 in G_ple:
            G_copy.add_edge(new_node, n2, data)
        elif n2 in nodes and n1 in G_cle:
            G_copy.add_edge(n1, new_node, data)