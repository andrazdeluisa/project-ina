import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from cdlib import algorithms
import infomap as imp
from collections import defaultdict
import cdlib as cdl
from copy import deepcopy
import random
from scipy.integrate import quad


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


def conductance(graph, partition):
    # doesn't work on graphs with removed nodes (indexing)
    partition = list(map(list, partition))
    sum_deg = [sum(list(map(lambda x: x[1], list(graph.degree(comm))))) for comm in partition[:10]]
    cuts = [nx.cut_size(graph, comm) for comm in partition[:10]]
    return [cuts[i] / sum_deg[i] for i in range(len(partition[:10]))]


github = nx.read_adjlist('networks/musae_git_edges.adj', create_using=nx.Graph, nodetype=int, delimiter=',')
random.seed(123)
print(len(max(nx.connected_components(github), key=len)) / 37700)
github_sub = github.subgraph(random.sample(list(github.nodes), 1000))
print(len(max(nx.connected_components(github_sub), key=len)) / 1000)

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
#comm_leid = algorithms.leiden(github)
#cond_leid = comm_leid.conductance(summary=False)
#print(comm_leid.conductance())

#comm_imp = infomap(github)
#cond_imp = comm_imp.conductance(summary=False)
#print(comm_imp.conductance())


# how to implement the removal? two options:
# a) remove top nodes, run community detection, evaluate the partition (e.g. conductance) -> other communities can form, track only the average conductance
# b) remove top nodes, keep same communities, evaluate partition -> we can track the dismantling of every single community
# in both we can choose to remove the nodes one by one or simultaneously

remove_pct = [0.001, 0.01, 0.05, 0.1]


def remove_and_detect(graph, centralities, community, remove_pct):
    N = graph.number_of_nodes()
    conductances = []
    avg_odfs = []
    cuts = []
    giant_cc = []
    for centrality in centralities:
        graph_tmp = nx.Graph(graph)
        rank = ranking(centrality(graph_tmp))
        partition = community(graph_tmp, epsilon=0.25)
        cond = [np.mean(partition.conductance(summary=False))]
        avg_odf = [np.mean(partition.avg_odf(summary=False))]
        cut = [np.mean(partition.cut_ratio(summary=False))]
        giant = [len(max(nx.connected_components(graph_tmp), key=len)) / N]
        for pct in remove_pct:
            nr_rem = int(pct * N)
            nodes = rank[:nr_rem]
            # warning: better to remove nodes from a copy of the graph
            graph_tmp.remove_nodes_from(nodes)
            partition_pct = community(graph_tmp, epsilon=0.25)
            cond_pct = np.mean(partition_pct.conductance(summary=False))
            cond.append(cond_pct)
            avg_odf.append(np.mean(partition_pct.avg_odf(summary=False)))
            cut.append(np.mean(partition_pct.cut_ratio(summary=False)))
            giant.append(len(max(nx.connected_components(graph), key=len)) / (N - nr_rem))
        conductances.append(cond)
        avg_odfs.append(avg_odf)
        cuts.append(cut)
        giant_cc.append(giant)
    return conductances, avg_odfs, cuts, giant_cc


def remove_and_keep(graph, centralities, community, remove_pct):
    N = graph.number_of_nodes()
    conductances = []
    avg_odfs = []
    cuts = []
    giant_cc = []
    for centrality in centralities:
        graph_tmp = nx.Graph(graph)
        rank = ranking(centrality(graph_tmp))
        partition = community(graph_tmp, epsilon=0.25)
        cond = [np.mean(partition.conductance(summary=False))]
        avg_odf = [np.mean(partition.avg_odf(summary=False))]
        cut = [np.mean(partition.cut_ratio(summary=False))]
        giant = [len(max(nx.connected_components(graph_tmp), key=len)) / N]
        for pct in remove_pct:
            nr_rem = int(pct * N)
            nodes = rank[:nr_rem]
            # warning: better to remove nodes from a copy of the graph
            graph_tmp.remove_nodes_from(nodes)
            cond_pct = np.mean(partition.conductance(summary=False))
            cond.append(cond_pct)
            avg_odf.append(np.mean(partition.avg_odf(summary=False)))
            cut.append(np.mean(partition.cut_ratio(summary=False)))
            giant.append(len(max(nx.connected_components(graph_tmp), key=len)) / (N - nr_rem))
        conductances.append(cond)
        avg_odfs.append(avg_odf)
        cuts.append(cut)
        giant_cc.append(giant)
    return conductances, avg_odfs, cuts, giant_cc


#print(remove_and_keep(github, [nx.degree_centrality, nx.eigenvector_centrality], algorithms.demon, remove_pct))


def find_best_nodes(graph):
    partitions = [algorithms.demon(graph, epsilon=eps) for eps in [0.25, 0.5, 0.75]]
    node_importance_area = defaultdict()
    node_importance_abs = defaultdict()
    sizes, conductances = ncp(partitions)
    x = range(2, max(sizes) + 1)
    for node in graph.nodes:
        # remove node 
        # compute conductance
        # compute ncp curve
        # compute difference in ncp
        # define node importance (area between ncp or abs diff)
        graph_tmp = nx.Graph(graph)
        graph_tmp.remove_nodes_from([node])
        partitions_tmp = deepcopy(partitions)
        for part in partitions_tmp:
            part.graph = graph_tmp
        sizes_tmp, conductances_tmp = ncp(partitions_tmp)
        f = lambda x: np.interp(x, sizes_tmp, conductances_tmp) - np.interp(x, sizes, conductances)
        node_importance_area[node] = quad(f, x[0], x[-1], points=x)
        node_importance_abs[node] = np.sign(max(f(x)) + min(f(x))) * max(abs(f(x)))
    return node_importance_area, node_importance_abs


def ncp(partitions):
    sizes = []
    conductances = []
    for partition in partitions:
        conductance = partition.conductance(summary=False)
        size = list(map(len, partition.communities))
        sizes += size
        conductances += conductance
    sizes_ord = np.argsort(sizes)
    sizes = list(np.array(sizes)[sizes_ord])
    conductances = list(np.array(conductances)[sizes_ord])
    res = []
    i = 0
    while i < len(sizes):
        tmp = [conductances[i]]
        j = i + 1
        while j < len(sizes) and sizes[i] == sizes[j]:
            tmp.append(conductances[j])
            j += 1
        res.append(min(tmp))
        i = j
    return np.unique(sizes), res


#karate = nx.karate_club_graph()
#partition = algorithms.demon(karate, epsilon=0.9, min_com_size=1)
#sizes, cond = ncp(partition)
#plt.plot(sizes, cond)
#plt.show()
partition2 = [algorithms.demon(github_sub, epsilon=eps, min_com_size=3) for eps in [0.25, 0.5, 0.75]]
sizes2, cond2 = ncp(partition2)
plt.plot(sizes2, cond2)
plt.show()

areas, diffs = find_best_nodes(github_sub)
# unfreeze graph
# remove nodes
# compute conductances
github_sub2 = nx.Graph(github_sub)
github_sub2.remove_nodes_from(random.sample(list(github_sub2.nodes), 50))
partition3 = deepcopy(partition2)
for part in partition3:
    part.graph = github_sub2
sizes3, cond3 = ncp(partition3)
plt.plot(sizes2, cond2)
plt.plot(sizes3, cond3)
plt.show()

x = list(range(2, max(sizes3) + 1))
y = np.interp(x, sizes3, cond3)
y1 = np.interp(x, sizes2, cond2)
f = lambda x: np.interp(x, sizes3, cond3) - np.interp(x, sizes2, cond2)
area3 = quad(f, x[0], x[-1], points=x)
abs_diff = max(abs(f(range(2,15))))


print('birc')