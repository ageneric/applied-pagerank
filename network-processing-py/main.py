import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt
from math import sqrt


DATA_EXPRESSIONS = '../network-data/expressions.csv'
DATA_TRRUST = '../network-data/TRRUST.csv'
DATA_STRING = '../network-data/STRING_by_gene.csv'

with open(DATA_EXPRESSIONS, 'r') as f_expressions:
    expressions = pd.read_csv(f_expressions)

with open(DATA_TRRUST, 'r') as f_trrust:
    trrust = pd.read_csv(f_trrust)

with open(DATA_STRING, 'r') as f_string_gene:
    string_gene = pd.read_csv(f_string_gene)


def get_gene_differential_expressions():
    gene_deg_totals = {}  # differentially expressed genes

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

        for gene in genes:
            if gene in gene_deg_totals:
                # we want to compute the mean of all gene
                gene_deg_totals[gene][0] += 1
                gene_deg_totals[gene][1] += row['rmeans_texpr'] - row['rmeans_nexpr']
            else:
                gene_deg_totals[gene] = [1, row['rmeans_texpr'] - row['rmeans_nexpr']]

    # compute the difference of the averages by dividing by the count
    ret = {}
    for key, value in gene_deg_totals.items():
        count, expression_difference_total = value
        ret[key] = abs(expression_difference_total / count)

    return ret


def get_gene_directed_neighbours(gene_name, checkdict):
    return trrust[(trrust['V1'] == gene_name) & (trrust['V2'].isin(checkdict))]


def GM(gene_a, gene_b, gene_deg):
    return sqrt(gene_deg[gene_a] * gene_deg[gene_b])


def RMS(gene_a, gene_b, gene_deg):
    return sqrt((gene_deg[gene_a]**2 + gene_deg[gene_b]**2) / 2)


def STRING(gene_a, gene_b, gene_deg):
    condition = (string_gene['gene1'] == gene_a) & (string_gene['gene2'] == gene_b)
    result = string_gene.loc[condition, 'combined_score']
    if not result.empty:
        return result.iloc[0] / 100  # n.b. result.iloc[0] gives you STRING score, not sure how to deal with it to make suitable weighting
    else:
        return RMS(gene_a, gene_b) # n.b. arbitrary rating for case where no STRING match is found


if __name__ == '__main__':
    network = nx.DiGraph()
    gene_deg = get_gene_differential_expressions()

    threshold = 2
    filtered_list = [x for x in gene_deg if abs(gene_deg[x]) > threshold]

    edge_weighting = STRING

    for gene in filtered_list:
        neighbours = get_gene_directed_neighbours(gene, filtered_list)
        neighbour_genes, neighbour_rel = neighbours['V2'], neighbours['V3']
        network.add_weighted_edges_from((n, gene, edge_weighting(n, gene, gene_deg)) for n in neighbour_genes)

    nx.draw(network, with_labels=True)
    plt.show()
