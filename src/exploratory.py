import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from cdlib import algorithms
import infomap as imp
from collections import defaultdict
import cdlib as cdl


def deg_distr(graph):
    degrees = np.asarray(list(dict(graph.degree).values()))
    _, pk = np.unique(degrees, return_counts=True)
    return pk


def ranking(centrality):
    ranked = sorted(centrality, key=centrality.get, reverse=True)
    return ranked


def infomap(graph):
    coms_to_node = defaultdict(list)

    im = imp.Infomap('')
    network = im.network()
    for e in graph.edges():
        network.addLink(int(e[0]), int(e[1]))
    im.run()
    for node in im.iterTree():
        if node.isLeaf():
            nid = node.physicalId
            module = node.moduleIndex()
            coms_to_node[module].append(nid)
   
    coms_infomap = [list(c) for c in coms_to_node.values()]
    return cdl.NodeClustering(coms_infomap, graph, "Infomap", method_parameters={"":""})


github = nx.read_adjlist('networks/musae_git_edges.adj', create_using=nx.Graph, nodetype=int, delimiter=',')

'''
print(nx.info(github))

pk = deg_distr(github)

plt.scatter(list(range(len(pk))), pk, s=2, c='b')
plt.yscale('log')
plt.xscale('log')
plt.title('Degree distribution')
plt.xlabel('k')
plt.ylabel('frequency')
plt.show()

# some centrality measures, add others
eig = nx.eigenvector_centrality(github)
eigrank = ranking(eig)
pr = nx.pagerank(github)
prrank = ranking(pr)
deg = nx.degree_centrality(github)
degrank = ranking(deg)
'''

# slower algorithms (check time complexity and estimate feasability)
# betw = nx.betweenness_centrality(github)
# betwrank = ranking(betw)
# clos = nx.closeness_centrality(github)
# closrank = ranking(clos)

# some community detection, add other algorithms
comm_leid = algorithms.leiden(github)
cond_leid = comm_leid.conductance(summary=False)
print(comm_leid.conductance())

comm_imp = infomap(github)
cond_imp = comm_imp.conductance(summary=False)
print(comm_imp.conductance())


# how to implement the removal? two options:
# a) remove top nodes, run community detection, evaluate the partition (e.g. conductance) -> other communities can form, track only the average conductance
# b) remove top nodes, keep same communities, evaluate partition -> we can track the dismantling of every single community
# in both we can choose to remove the nodes one by one or simultaneously

remove_pct = [0.01, 0.05, 0.1, 0.2]


def remove(graph, centralities, community, remove_pct):
    N = graph.number_of_nodes()
    conductances = []
    for centrality in centralities:
        rank = ranking(centrality(graph))
        cond = [np.mean(community(graph).conductance(summary=False))]
        for pct in remove_pct:
            nr_rem = int(pct * N)
            nodes = rank[:nr_rem]
            # warning: better to remove nodes from a copy of the graph
            graph.remove_nodes_from(nodes)
            cond_pct = np.mean(community(graph).conductance(summary=False))
            cond.append(cond_pct)
        conductances.append(cond)
    return conductances


print(remove(github, [nx.degree_centrality, nx.pagerank], infomap, remove_pct))
