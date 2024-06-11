import networkx as nx
import numpy as np


def draw_network_pagerank(network, top_gene_names, top_pageranks, top_n=40):
    view_network = nx.subgraph(network, top_gene_names[:top_n])
    scale = 900 / max(0.01, top_pageranks[0])
    nx.draw(view_network, node_size=np.array(top_pageranks[:top_n]) * scale, with_labels=True)

