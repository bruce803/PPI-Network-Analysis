import networkx as nx
import pylab as plt
import extract
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

[n_layer0,n_layer1,n_layer2,n_layer3,edge_all] = extract.extract_node_edges_from_file(cle_up, ple_up, g, nodes_inv_dict, 3)
node_all = sum([n_layer0, n_layer1, n_layer2,n_layer3], [])



# plottting
g_pleuucle_1hop = nx.Graph()
g_pleuucle_1hop.add_nodes_from(node_all)
g_pleuucle_1hop.add_edges_from(edge_all)

lens = [len(n_layer0), len(n_layer1), len(n_layer2)]

pos_0 = {}
x0 = 1
const = 2
y = 100
for i in range(len(n_layer0)):
    pos_0[n_layer0[i]] = [x0, y - i * (max(lens) / len(n_layer0))]

x1 = 3
# y2=100
pos_1 = {}
for j in range(len(n_layer1)):
    pos_1[n_layer1[j]] = [x1, y - j * (max(lens) / len(n_layer1))]

x2 = 7
# y2=100
pos_2 = {}
for n in range(len(n_layer2)):
    pos_2[n_layer2[n]] = [x2, y - n * (max(lens) / len(n_layer2))]

x3=11
# y3=700
pos_3={}
for m in range(len(n_layer3)):
    pos_3[n_layer3[m]] = [x3, y - m * (max(lens) / len(n_layer3))]

# nx.draw(g_pleddcle_1hop, pos=position_MultiPartiteGraph(g_pleddcle_1hop, ['start', 'mid', 'end']))

nx.draw_networkx_nodes(g_pleuucle_1hop, pos_0, nodelist=n_layer0, node_color='g', node_size=200, alpha=0.8)
nx.draw_networkx_nodes(g_pleuucle_1hop, pos_1, nodelist=n_layer1, node_color='r', node_size=200, alpha=0.8)
nx.draw_networkx_nodes(g_pleuucle_1hop, pos_2, nodelist=n_layer2, node_color='r', node_size=100, alpha=0.8)
nx.draw_networkx_nodes(g_pleuucle_1hop, pos_3, nodelist=n_layer3, node_color='c', node_size=50, alpha=0.8)

# edges
pos = {}
pos.update(pos_0)
pos.update(pos_1)
pos.update(pos_2)
pos.update(pos_3)
nx.draw_networkx_edges(g_pleuucle_1hop, pos, edgelist=nx.edges(g_pleuucle_1hop), width=1, alpha=0.8, edge_color='k')
nx.draw_networkx_labels(g_pleuucle_1hop, pos, font_size=8, font_family='sans-serif')

fig = plt.gcf()
plt.show()
plt.draw()
fig.set_size_inches(8, 5)
plt.axis('off')
# plt.title('The 2-hops subnetwork from CLE up to PLE up')
# fig.savefig('g_pleuucle_2hop_new.png', dpi=900)

