import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from cdlib import algorithms
import infomap as imp
from collections import defaultdict
import cdlib as cdl
from copy import deepcopy
import random
from time import time


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
random.seed(123)


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


def find_best_nodes(graph, top_nodes):
    partitions = [algorithms.demon(graph, epsilon=eps) for eps in [0.25]]#, 0.5, 0.75]]
    node_importance_area = defaultdict()
    node_importance_abs = defaultdict()
    sizes = []
    conductances = []
    edges_insides = []
    communities = []
    for partition in partitions:
        sizes += partition.size(summary=False)
        conductances += partition.conductance(summary=False)
        edges_insides += partition.edges_inside(summary=False)
        communities += partition.communities
    sizes_ncp, conductances_ncp = ncp(sizes, conductances)
    with open('src/github_ncp.csv', 'w+') as file:
        file.write('Size,Cond\n')
        for i in range(len(sizes_ncp)):
            file.write(','.join([str(sizes_ncp[i]), str(conductances_ncp[i]) + '\n']))
        file.close()
    plt.plot(sizes_ncp, conductances_ncp)
    plt.savefig('report/github.png', dpi=600)
    plt.close()
    x = range(2, max(sizes_ncp) + 1)
    values = np.interp(x, sizes_ncp, conductances_ncp)
    print('partitions computed')
    i = 0
    for node in top_nodes:
        # remove node 
        # compute conductance
        # compute ncp curve
        # define node importance (area between ncp or abs diff)
        sizes_tmp, conductances_tmp = conductance_fast(graph, communities, conductances, edges_insides, node)
        sizes_tmp, conductances_tmp = ncp(sizes_tmp, conductances_tmp)
        f_values = np.interp(x, sizes_tmp, conductances_tmp) - values
        node_importance_area[node] = integrate(f_values)
        node_importance_abs[node] = np.sign(max(f_values) + min(f_values)) * max(abs(f_values))
        i += 1
        if i % 10 == 0:
            print(i)
    return dict(node_importance_area), dict(node_importance_abs)


def ncp(sizes, conductances):
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
    return list(np.unique(sizes)), res


def integrate(y):
    area = 0
    for i in range(len(y) - 1):
        if y[i] * y[i+1] == 0:
            area += y[i] / 2 + y[i+1] / 2
        elif np.sign(y[i]) == np.sign(y[i+1]):
            area += np.sign(y[i]) * (min(abs(y[i]), abs(y[i+1])) + abs(y[i+1] - y[i]) / 2)
        else:
            t = - y[i] / (y[i+1] - y[i])
            area += (t * y[i] + (1 - t) * y[i+1]) / 2
    return area


def conductance_fast(graph, communities, conductances, edges_inside, node):
    sizes = list(map(len, communities))
    contained = [node in comm for comm in communities]
    neighb = list(graph.neighbors(node))
    contained_neighb = []
    edges_outside = []
    corr_conductances = [0 for _ in range(len(conductances))]
    for j in neighb:
        contained_neighb.append([j in comm for comm in communities])
    for i in range(len(conductances)):
        edges_outside.append(sum(np.array(graph.degree)[communities[i], 1]))
        edges_outside[i] -= 2 * edges_inside[i]
        n_inside = 0
        for tmp in contained_neighb:
            if tmp[i]:
                n_inside += 1
        if contained[i]:
            sizes[i] -= 1
            n_outside = len(neighb) - n_inside
        else:
            n_outside = n_inside
            n_inside = 0
        corr_conductances[i] = (conductances[i] * (2 * edges_inside[i] + edges_outside[i]) - n_outside) / (2 * (edges_inside[i] - n_inside) + edges_outside[i] - n_outside)
    return sizes, corr_conductances


#karate = nx.karate_club_graph()
#dolphins = nx.read_adjlist('networks/dolphins.adj', create_using=nx.Graph, nodetype=int)

eig = nx.eigenvector_centrality(github)
rank = ranking(eig)

start = time()
areas, diffs = find_best_nodes(github, rank[:100])
print(time() - start)

with open('src/tmp.csv', 'w+') as file:
    file.write('Id,Area,Diff\n')
    for i in rank[:100]:
        file.write(','.join([str(i), str(areas[i]), str(diffs[i])]) + '\n')
    file.close()

# if graph is full (contains all nodes from 1 to n)
areas = list(areas.values())
diffs = list(diffs.values())
# otherwise take care of indexing in dict (built-in centrality measures also return dict)

print('birc')