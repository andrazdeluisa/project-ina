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
