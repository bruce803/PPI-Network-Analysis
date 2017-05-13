import itertools
import argparse
import os



# load gene list
def readGeneFromFile(file_path):
    gene_list = []
    with open(file_path, "r") as fp:
        for line in fp:
            tmp = line.rstrip('\n')
            try:
                gene_list.append(tmp)
            except:pass
    return gene_list


def geneListInter(list1, list2):
    '''get the instersection of two sets of nodes or edges. The nodes (char) or edges(tuple) should be stored as List.
    For example list1=geneList, and list2=netNodes are genes' names'''
    return list(set(list1)&set(list2))


# convert from NAME to id number ('key' in nodes_dict)
def name2id(inter_immune,nodes_inv_dict):
    inter_immune_id = []
    for i in range(len(inter_immune)):
        inter_immune_id.append(nodes_inv_dict[inter_immune[i]])
    return inter_immune_id


def complete_bipartite(list1,list2):
    '''generate a complete bipartite graph between list1 and list2'''
    c12 = list(itertools.product(list1, list2)) #考虑(1,4)和（4，1）两种情况
    c21 = list(itertools.product(list2, list1))
    return c12 + c21


def read_dict_from_folder(indir):
    '''read the function or pathway into dictionary, key:pathway name, value:pathway'''
    pathFunDict = {}
    for root, dirs, filenames in os.walk(indir):
        for fn in filenames:
            path = os.path.join(indir,fn)
            geneList = readGeneFromFile(path)
            pathFunDict[fn]= geneList
    return pathFunDict


def str2bool(v):
    '''parse boolen values with argparse'''
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
