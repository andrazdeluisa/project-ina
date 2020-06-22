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


github = nx.read_adjlist('networks/musae_git_edges.adj', create_using=nx.Graph, nodetype=int, delimiter=',')
random.seed(123)
# print(nx.info(github))
# pk = deg_distr(github)
#
#plt.scatter(list(range(len(pk))), pk, s=2, c='b')
#plt.yscale('log')
#plt.xscale('log')
#plt.title('Degree distribution')
#plt.xlabel('k')
#plt.ylabel('frequency')
#plt.show()


def find_best_nodes(graph, top_nodes):
    partitions = [algorithms.demon(graph, epsilon=eps) for eps in [0.25, 0.5]]#, 0.75]]
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
areas, diffs = find_best_nodes(github, rank[:200])
print(time() - start)

with open('src/tmp.csv', 'w+') as file:
    file.write('Id,Area,Diff\n')
    for i in rank[:200]:
        file.write(','.join([str(i), str(areas[i]), str(diffs[i])]) + '\n')
    file.close()

print('birc')