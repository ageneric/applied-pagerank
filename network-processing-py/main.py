import pandas
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt

from gene_data import get_gene_differential_expressions
from weighting import WeightVectorMethod, expressions, get_TFLink_subset, get_STRING_subset, TF, TARGET
from pagerank import linear_pagerank
from graphic import draw_network_pagerank


# filter for the top differentially expressed genes
# 6.9 = keep top ~90%. 14.6 = keep top ~80%.
# don't recommend dropping more than 20% as TP53 (important) has low expression
THRESHOLD_DEG = 14.6

def write_to_graphml():
    import pathlib
    nx.write_graphml(network, pathlib.Path.cwd() / 'graph.graphml')


if __name__ == '__main__':
    print('Imported modules.')

    gene_deg = get_gene_differential_expressions(expressions)
    filter_list = [g for (g, differential) in gene_deg.items()
                   if abs(differential) > THRESHOLD_DEG]
    tflink_df = get_TFLink_subset(filter_list)
    string_df = get_STRING_subset(filter_list)
    print('Imported datasets.')

    df = tflink_df.merge(string_df, on=(TARGET, TF))
    weighting = WeightVectorMethod(gene_deg, df).STRING
    print('Merged datasets.')

    print(f'''Genes with known expressions               {len(gene_deg)}
Genes passing threshold expression         {len(filter_list)}
 -- and their interactions in TFLink       {len(tflink_df)}
 -- and their interactions in STRING       {len(string_df)}
Interactions in TFLink + STRING merged     {len(df)}
Weighting method                           {weighting.__name__}''')

    # Free up memory from the separate datasets
    del tflink_df
    del string_df
    network = nx.DiGraph()

    # Since few genes have out-links, we filter down to only the genes with
    # out-links and iterate through them. This is faster than enumerating
    # through all genes by a significant margin.
    genes_with_out_links = df[TF].unique()

    for i, gene in enumerate(genes_with_out_links):
        neighbours = df[df[TF] == gene]
        if not neighbours.empty:
            weights = weighting(gene, neighbours)
            print(i, end=' ')
            network.add_weighted_edges_from((neighbour, gene, weight) for neighbour, weight in zip(neighbours[TARGET], weights))
            print(gene)
        if i % 1000 == 999:
            print(f'Generating NetworkX graph: {i} complete')
    print('Generated NetworkX graph.')

    stochastic_network = nx.stochastic_graph(network)

    matrix_P = nx.adjacency_matrix(stochastic_network, stochastic_network.nodes,
                                   weight='weight').todense()

    pagerank = linear_pagerank(matrix_P)
    pagerank_dict = {k: v for k, v in zip(network.nodes, pagerank)}
    nx.set_node_attributes(network, pagerank_dict, name='pagerank')

    genes = [(item[0], item[1]) for item in zip(network.nodes, pagerank)]
    top_genes = sorted(genes, key=lambda x: -x[1])
    top_gene_names = [data[0] for data in top_genes]
    top_pageranks = [data[1] for data in top_genes]
    print('Computed PageRank list.')

    # Visualise 'top 50' genes
    draw_network_pagerank(network, top_gene_names, top_pageranks)
    plt.show()

