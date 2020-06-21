import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from cdlib import algorithms
import infomap as imp
from collections import defaultdict
import cdlib as cdl
from copy import deepcopy
from networkx.generators.community import LFR_benchmark_graph
from exploratory import infomap


benchmark = LFR_benchmark_graph(n=20, tau1=2.1, tau2=1.5, mu=0.7, average_degree=5, seed=1)
er = nx.gnm_random_graph(n=100, m=500)
aaa = infomap(er)

communities = {frozenset(benchmark.nodes[v]['community']) for v in benchmark}

print(len(communities))
print(communities)

tmp = nx.draw_networkx(benchmark)
plt.show(tmp)

tmp = np.array(benchmark.degree)
communities_list = list(map(list, communities))
sum_deg = list(map(sum, [tmp[comm, 1] for comm in communities_list]))
cuts = [nx.cut_size(benchmark, comm) for comm in communities_list]
conductance = [cuts[i] / sum_deg[i] for i in range(len(communities))]


def conductance(graph, partition):
    partition = list(map(list, partition))
    sum_deg = list(map(sum, [np.array(graph.degree)[comm, 1] for comm in partition]))
    cuts = [nx.cut_size(graph, comm) for comm in partition]
    return[cuts[i] / sum_deg[i] for i in range(len(partition))]


nx.community.quality.performance(benchmark, communities)

print('birc')
