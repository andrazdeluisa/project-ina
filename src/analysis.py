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
import pandas as pd


def ranking(centrality):
    ranked = sorted(centrality, key=centrality.get, reverse=True)
    return ranked


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


#github_ncp = pd.read_csv('src/github_ncp.csv', sep=',')

#plt.plot(github_ncp['Size'], github_ncp['Cond'])
#plt.xscale('log')
#plt.xlabel('Community size')
#plt.ylabel('Conductance')
#plt.title('Estimation of NPC plot in github network')
#plt.savefig('report/github_npc.png', dpi=600)

# Degree centrality measures

github = nx.read_adjlist('networks/musae_git_edges.adj', create_using=nx.Graph, nodetype=int, delimiter=',')

eig = nx.eigenvector_centrality(github)
eigrank = ranking(eig)
#pr = nx.pagerank(github)
#prrank = ranking(pr)
deg = nx.degree_centrality(github)
degrank = ranking(deg)

github200 = np.array(pd.read_csv('src/github200.csv', sep=','))

#area = list(np.array(github200['Area'])[np.argsort(github200['Area'])[::-1]])
arearank = github200[np.argsort(github200[:,1])[::-1], 0]
#abs = list(np.array(github200['Diff'])[np.argsort(github200['Diff'])[::-1]])
absrank = github200[np.argsort(github200[:,2])[::-1], 0]
abs_list = [github200[github200[:, 0] == node, 2][0] for node in eigrank[:200]]
area = [github200[github200[:, 0] == node, 1][0] for node in eigrank[:200]]
eig_list = [eig[node] for node in eigrank[:200]]
deg_list = [deg[node] for node in eigrank[:200]]
print(np.corrcoef(abs_list, eig_list))
print(np.corrcoef(area, eig_list))
print(np.corrcoef(abs_list, area))
print(np.corrcoef(abs_list, deg_list))
print(np.corrcoef(area, deg_list))
print(np.corrcoef(deg_list, eig_list))

'''
partitions = [algorithms.demon(github, epsilon=eps) for eps in [0.25, 0.5]]
sizes = []
conductances = []
for partition in partitions:
    sizes += partition.size(summary=False)
    conductances += partition.conductance(summary=False)
sizes_ncp, conductances_ncp = ncp(sizes, conductances)
graph_tmp = nx.Graph(github)
graph_tmp.remove_nodes_from(arearank[:20])
sizes2 = []
conductances2 = []
for partition in partitions:
    partition.graph = graph_tmp
    sizes2 += partition.size(summary=False)
    conductances2 += partition.conductance(summary=False)
sizes2_ncp, conductances2_ncp = ncp(sizes2, conductances2)
with open('src/github2_ncp.csv', 'w+') as file:
    file.write('Size,Cond\n')
    for i in range(len(sizes_ncp)):
        file.write(','.join([str(sizes_ncp[i]), str(conductances_ncp[i]) + '\n']))
    for j in range(len(sizes2_ncp)):
        file.write(','.join([str(sizes2_ncp[j]), str(conductances2_ncp[j]) + '\n']))
    file.close()
plt.plot(sizes_ncp, conductances_ncp, c='b')
plt.plot(sizes2_ncp, conductances2_ncp, c='r')
plt.xscale('log')
plt.xlabel('Community size')
plt.ylabel('Conductance')
plt.title('NCP plot before and after removing top 20 nodes')
plt.legend(['Original network', 'Removed nodes'])
plt.savefig('report/github2.png', dpi=600)
plt.close()
'''

print('birc')