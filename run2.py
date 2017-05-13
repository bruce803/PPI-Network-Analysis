import networkx as nx
import pylab as plt
import argparse
import merge
import extract
import utils
try:
    import cPickle as pickle
except ImportError:
    import pickle

#load network
g = nx.read_gpickle('graph.gpickle')

# load nodes dictionary
with open('nodeDict.data','rb') as f:
    nodes_dict = pickle.load(f)

# load inverse nodes dictionary
with open('nodeDictInv.data','rb') as f:
    nodes_inv_dict = pickle.load(f)

# load gene list
cle_up = 'C:\\Users\\mrliu\\Documents\\copd-net\\data\\paren_CLE_C_Lmlme_result_up_symbol.txt'
ple_up = 'C:\\Users\\mrliu\\Documents\\copd-net\\data\\paren_PLE_C_Lmlme_result_up_symbol.txt'
pathwayDir = 'C:\\Users\\mrliu\\Documents\\copd-net\\data\\pathway_list'
funDir = 'C:\\Users\\mrliu\\Documents\\copd-net\\data\\function_list'

funDict = utils.read_dict_from_folder(funDir) #396
pathwayDict = utils.read_dict_from_folder(pathwayDir) #1807


def experiment1(gene1ist1, gene1ist2, graph, nodesinvdict):
    [n_layer0, n_layer1, edge_all] = extract.extract_node_edges_from_file(gene1ist1, gene1ist2, graph, nodesinvdict, 1)
    # node_all = sum([n_layer0, n_layer1], [])
    return n_layer0, n_layer1, edge_all


def experiment2(gene1ist1, gene1ist2, graph, nodesinvdict):
    [n_layer0, n_layer1, n_layer2, edge_all] = extract.extract_node_edges_from_file(gene1ist1, gene1ist2, graph, nodesinvdict, 2)
    node_all = sum([n_layer0, n_layer1, n_layer2], [])

    g_2hop= nx.Graph()
    g_2hop.add_nodes_from(node_all)
    g_2hop.add_edges_from(edge_all)

    node_02_copy = n_layer0 + n_layer2  # 37
    g_2hop_copy = nx.Graph()
    g_2hop_copy.add_nodes_from(node_02_copy)


    inter_node = merge.node_inter_g(pathwayDict, n_layer1) #613, pathway
    # inter_node = merge.node_inter_g(funDict, n_layer1) #67, function

    for k, v in inter_node.items():
        #     merge_id = [nodes_inv_dict[gene_name] for gene_name in v] # display with id
        merge.copy_merge_nodes(g_2hop, g_2hop_copy, n_layer0, n_layer2, v, k)



    node_layer1_copy = [item for item in g_2hop_copy.nodes() if item not in node_02_copy] # g_3hop_copy.nodes - node_copy

    lens = [len(n_layer0), len(node_layer1_copy), len(n_layer2)]

    pos_0 = {}
    x0 = 1
    const = 2
    y = 100
    for i in range(len(n_layer0)):
        pos_0[n_layer0[i]] = [x0, y - i * (max(lens) / len(n_layer0))]

    x1 = 3
    pos_1 = {}
    for j in range(len(node_layer1_copy)):
        pos_1[node_layer1_copy[j]] = [x1, y - j * (max(lens) / len(node_layer1_copy))]

    x2 = 5
    pos_2 = {}
    for n in range(len(n_layer2)):
        pos_2[n_layer2[n]] = [x2, y - n * (max(lens) / len(n_layer2))]

    nx.draw_networkx_nodes(g_2hop_copy, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
    nx.draw_networkx_nodes(g_2hop_copy, pos_1, nodelist=node_layer1_copy, node_color='r', node_size=300, alpha=0.8)
    nx.draw_networkx_nodes(g_2hop_copy, pos_2, nodelist=n_layer2, node_color='c', node_size=300, alpha=0.8)
    # nx.draw_networkx_nodes(g_pleuucle_1hop,pos_3,nodelist=n_layer3,node_color='c',node_size=300,alpha=0.8)

    # edges
    pos = {}
    pos.update(pos_0)
    pos.update(pos_1)
    pos.update(pos_2)
    # pos.update(pos_3)
    nx.draw_networkx_edges(g_2hop_copy, pos, edgelist=nx.edges(g_2hop_copy), width=1, alpha=0.8,
                           edge_color='k')
    nx.draw_networkx_labels(g_2hop_copy, pos, font_size=8, font_family='sans-serif')

    fig = plt.gcf()
    plt.show()
    plt.draw()
    fig.set_size_inches(8, 5)
    plt.axis('off')
    plt.title('The merged subnetwork from CLE up to PLE up')
    # fig.savefig('g_pathway_merge.png', dpi=900)
    # fig.savefig('g_function_merge.png', dpi=900)


def experiment3(gene1ist1, gene1ist2, graph, nodesinvdict, mergeBoth):
    [n_layer0, n_layer1, n_layer2, n_layer3, edge_all] = extract.extract_node_edges_from_file(gene1ist1, gene1ist2, graph,
                                                                                              nodesinvdict, 3)
    node_all = sum([n_layer0, n_layer1, n_layer2, n_layer3], [])

    # return n_layer0, n_layer1, n_layer2, n_layer3, edge_all

    g_3hop= nx.Graph()
    g_3hop.add_nodes_from(node_all)
    g_3hop.add_edges_from(edge_all)

    node_03_copy = n_layer0 + n_layer3  # 37
    g_3hop_copy = nx.Graph()
    g_3hop_copy.add_nodes_from(node_03_copy)

    # return n_layer0, n_layer1, n_layer2, n_layer3, edge_all, g_3hop_copy,

    if mergeBoth:
        # True
        inter_node = merge.node_both_inter_g(pathwayDict, n_layer1, n_layer2) #198, pathway
        # inter_node = merge.node_both_inter_g(funDict, n_layer1, n_layer2)  # 28, function
    elif not mergeBoth:
        n_layer12 = n_layer1 + n_layer2
        # #False
        inter_node = merge.node_inter_g(pathwayDict, n_layer12) #613, pathway
        # inter_node = merge.node_inter_g(funDict, n_layer12) #67, function


    for k, v in inter_node.items():
        #     merge_id = [nodes_inv_dict[gene_name] for gene_name in v] # display with id
        merge.copy_merge_nodes(g_3hop, g_3hop_copy, n_layer0, n_layer3, v, k)



    node_layer1_copy = [item for item in g_3hop_copy.nodes() if item not in node_03_copy] # g_3hop_copy.nodes - node_copy

    lens = [len(n_layer0), len(node_layer1_copy), len(n_layer3)]

    pos_0 = {}
    x0 = 1
    const = 2
    y = 100
    for i in range(len(n_layer0)):
        pos_0[n_layer0[i]] = [x0, y - i * (max(lens) / len(n_layer0))]

    x1 = 3
    pos_1 = {}
    for j in range(len(node_layer1_copy)):
        pos_1[node_layer1_copy[j]] = [x1, y - j * (max(lens) / len(node_layer1_copy))]

    x2 = 5
    pos_2 = {}
    for n in range(len(n_layer3)):
        pos_2[n_layer3[n]] = [x2, y - n * (max(lens) / len(n_layer3))]

    nx.draw_networkx_nodes(g_3hop_copy, pos_0, nodelist=n_layer0, node_color='g', node_size=300, alpha=0.8)
    nx.draw_networkx_nodes(g_3hop_copy, pos_1, nodelist=node_layer1_copy, node_color='r', node_size=300, alpha=0.8)
    nx.draw_networkx_nodes(g_3hop_copy, pos_2, nodelist=n_layer3, node_color='c', node_size=300, alpha=0.8)
    # nx.draw_networkx_nodes(g_pleuucle_1hop,pos_3,nodelist=n_layer3,node_color='c',node_size=300,alpha=0.8)

    # edges
    pos = {}
    pos.update(pos_0)
    pos.update(pos_1)
    pos.update(pos_2)
    # pos.update(pos_3)
    nx.draw_networkx_edges(g_3hop_copy, pos, edgelist=nx.edges(g_3hop_copy), width=1, alpha=0.8,
                           edge_color='k')
    nx.draw_networkx_labels(g_3hop_copy, pos, font_size=8, font_family='sans-serif')

    fig = plt.gcf()
    plt.show()
    plt.draw()
    fig.set_size_inches(8, 5)
    plt.axis('off')
    plt.title('The merged subnetwork from CLE up to PLE up')
    # fig.savefig('g_pathway_merge.png', dpi=900)
    # fig.savefig('g_function_merge.png', dpi=900)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PPI Network Analysis')

    parser.add_argument('--layer', '-l', type=int, choices=[1, 2, 3], default=3)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--merge", "-m", action="store_true", default=False)
    group.add_argument('--delete', '-d', nargs='+', type=int)
    parser.add_argument('--pvalue', '-p', action="store", dest="fdr", type=float)
    args = parser.parse_args()

    if args.merge:
        if args.layer ==1:
            experiment1(cle_up, ple_up, g, nodes_inv_dict)
        elif args.layer == 2:
            experiment2(cle_up, ple_up, g, nodes_inv_dict)
        elif args.layer == 3:
            experiment3(cle_up, ple_up, g, nodes_inv_dict, mergeBoth=True)
            print('merge the middle nodes with fdr {0}'.format(args.fdr))
    elif not args.merge:
        if args.layer ==1:
            experiment1(cle_up, ple_up, g, nodes_inv_dict)
        elif args.layer == 2:
            experiment2(cle_up, ple_up, g, nodes_inv_dict)
        elif args.layer == 3:
            experiment3(cle_up, ple_up, g, nodes_inv_dict, mergeBoth=False)
            print('merge both the two layer node with fdr {0}'.format(args.fdr))

    # if args.delete:
    #     print('delte the node {0}'.format(args.delete))
