import pandas
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt

from gene_data import get_gene_differential_expressions
from weighting import WeightVectorMethod, expressions, get_TRRUST_subset, \
    get_TFLink_subset, get_STRING_subset, TF, TARGET
from pagerank import linear_system_pagerank, get_personalisation_vector_by_deg, \
    format_pagerank
from graphic import draw_network_pagerank


# filter for the top differentially expressed genes
# 6.9 = keep top ~90%. 14.6 = keep top ~80%.
# don't recommend dropping more than 20% as TP53 (important) has low expression
THRESHOLD_DEG = -1


def unzip(iterable):
    return zip(*iterable)

def generate_networkx_graph(_df, _weighting):
    network = nx.DiGraph()

    # Since few genes have out-links, we filter down to only the genes with
    # out-links and iterate through them. This is faster than enumerating
    # through all genes by a significant margin.
    genes_with_out_links = _df[TF].unique()
    num_genes_with_out_links = len(genes_with_out_links)

    for i, gene in enumerate(genes_with_out_links):
        neighbours = _df[_df[TF] == gene]
        if not neighbours.empty:
            weights = _weighting(gene, neighbours)
            network.add_weighted_edges_from((neighbour, gene, weight)
                                            for neighbour, weight in zip(neighbours[TARGET], weights))
        if i % 100 == 99:
            print(f'Generating NetworkX graph: {i/num_genes_with_out_links:.2%} complete')

    return network


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

    network = generate_networkx_graph(df, weighting)
    print('Generated NetworkX graph.')

    stochastic_network = nx.stochastic_graph(network)
    matrix_P = nx.adjacency_matrix(stochastic_network, stochastic_network.nodes,
                                   weight='weight').todense()
    personalisation = get_personalisation_vector_by_deg(network.nodes, gene_deg)
    print('Formulated matrix problem.')

    pagerank = linear_system_pagerank(matrix_P, v=personalisation, alpha=0.85)
    print('Computed PageRanks.')

    pagerank_dict, top_genes = format_pagerank(pagerank, network)
    top_gene_names, top_pageranks = unzip(top_genes)

    # Visualise 'top 50' genes
    draw_network_pagerank(network, pagerank_dict, top_gene_names, top_pageranks)
