import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt


DATA_EXPRESSIONS = '../network-data/expressions.csv'
DATA_TRRUST = '../network-data/TRRUST.csv'

with open(DATA_EXPRESSIONS, 'r') as f_expressions:
    expressions = pd.read_csv(f_expressions)

with open(DATA_TRRUST, 'r') as f_trrust:
    trrust = pd.read_csv(f_trrust)


def get_gene_differential_expressions():
    gene_deg = {}  # differentially expressed genes

    # group together and process the expression values by gene
    # iterrows() is slow, but we will only iterate over ~20000 samples
    for row in expressions.iterrows():
        row = row[1]
        if str(row['gene']) == 'nan':
            continue

        if '///' in row['gene']:
            genes = row['gene'].split(' /// ')
        else:
            genes = [row['gene']]

        differential = row['rmeans_texpr'] - row['rmeans_nexpr']

        for gene in genes:
            if gene in gene_deg:
                gene_deg[gene] += differential
            else:
                gene_deg[gene] = differential

    return gene_deg


def get_gene_directed_neighbours(gene_name):
    return trrust[trrust['V1'] == gene_name]


if __name__ == '__main__':
    network = nx.DiGraph()
    gene_deg = get_gene_differential_expressions()

    # TODO: generate graph directions and weights. This code draws the few
    # genes with greatest differential expression (and their connections)
    # and doesn't consider direction or weight

    for gene in gene_deg:
        if abs(gene_deg[gene]) > 5:
            neighbours = get_gene_directed_neighbours(gene)
            neighbour_genes, neighbour_rel = neighbours['V2'], neighbours['V3']
            network.add_edges_from((gene, n) for n in neighbour_genes)

    nx.draw(network, with_labels=True)
    plt.show()
