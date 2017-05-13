import networkx as nx
import itertools
try:
    import cPickle as pickle
except ImportError:
    import pickle

#load network

imm_net = 'C:\\Users\\mrliu\\Documents\\copd-net\\data\\immune_genes_human.txt'
#store all the nodes (ids,names) into Dict
nodes = []
with open(imm_net, "r") as fp:
    for line in itertools.islice(fp, 1, 1887): # read from row 2 to row 1887
        tmp = line.rstrip('\n') # delete change line
        tmp = tmp.split(" ")  # split with backspace
        try:
            nodes.append((int(tmp[0]),tmp[1].replace('\"','')))
        except:pass
nodes_dict = dict(nodes)


# serialize the nodes (id:name) dictionary
with open('nodeDict.data', 'wb') as f:
    pickle.dump(nodes_dict, f)

# store edges into tuple list
# num_lines = sum(1 for line in open(c)) # count the number of lines of network file 13773
edges = []
with open(imm_net, "r") as fp:
    for line in itertools.islice(fp, 1889, 13773): # read from row 1889 to end
        tmp = line.rstrip('\n 1')
        tmp = tmp.split(" ")
        try:
            edges.append((int(tmp[0]),int(tmp[1])))
        except:pass

# covert (4,33) to ('UBC', 'ITCH')
edges_id = [(nodes_dict[inNode], nodes_dict[outNode]) for inNode, outNode in edges]
#whole graph
#load nodes from Dictionary,g.size()=g.number_of_edges()=10821, g.number_of_nodes()=1886
node_id = nodes_dict.values() # the id are like this, 'MAPK8','MAP2K7', 'UBC', 'ITCH'.....
g = nx.Graph()
# g.add_nodes_from(node_id,name=nodes_dict.values())
g.add_nodes_from(node_id)
g.add_edges_from(edges_id)

# write the graph to file

nx.write_graphml(g,"graph.graphml")
nx.write_gml(g,"graph.gml")
nx.write_gpickle(g,"graph.gpickle")



# invert the nodes mapping
nodes_inv_dict = dict((v, k) for k, v in nodes_dict.items())
# nodes_inv_dict['XCL2'] #988


# serialize the inverse dictionary
with open('nodeDictInv.data', 'wb') as f:
    pickle.dump(nodes_inv_dict, f)