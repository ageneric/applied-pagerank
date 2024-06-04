import pandas as pd
from math import sqrt


DATA_EXPRESSIONS = '../network-data/expressions.csv'
DATA_TRRUST = '../network-data/TRRUST.csv'
DATA_STRING = '../network-data/STRING_by_gene.csv'

with open(DATA_EXPRESSIONS, 'r') as f_expressions:
    expressions = pd.read_csv(f_expressions)

with open(DATA_TRRUST, 'r') as f_trrust:
    trrust = pd.read_csv(f_trrust)


def get_STRING_subset(filter_list):
    with open(DATA_STRING, 'r') as f_string_gene:
        string_everything = pd.read_csv(f_string_gene)

    # Prepare STRING by only keeping edges between genes in the filter list
    return string_everything[(string_everything['name1'].isin(filter_list))
                             & (string_everything['name2'].isin(filter_list))]


class WeightMethod:
    DEFAULT_WEIGHT_NO_STRING = 0.1

    def __init__(self, deg, STRING_gene_data):
        self.deg = deg
        self.STRING_gene_data = STRING_gene_data

    def GM(self, gene_a, gene_b):
        return sqrt(self.deg[gene_a] * self.deg[gene_b])

    def RMS(self, gene_a, gene_b):
        return sqrt((self.deg[gene_a]**2 + self.deg[gene_b]**2) / 2)

    def STRING(self, gene_a, gene_b):
        condition = ((self.STRING_gene_data['gene1'] == gene_a)
                     & (self.STRING_gene_data['gene2'] == gene_b))
        result = self.STRING_gene_data.loc[condition, 'combined_score']
        if not result.empty:
            # n.b. result.iloc[0] gives you STRING score in thousandths
            return result.iloc[0] / 1000
        else:
            # n.b. arbitrary rating for case where no STRING match is found
            return self.DEFAULT_WEIGHT_NO_STRING

    def product_STRING(self, gene_a, gene_b):
        # Use STRING weight and multiply by RMS
        return self.RMS(gene_a, gene_b) * self.STRING(gene_a, gene_b)
