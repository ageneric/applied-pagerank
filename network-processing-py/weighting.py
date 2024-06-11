import pandas as pd
import numpy as np
from math import sqrt


DATA_EXPRESSIONS = '../network-data/expressions.csv'
DATA_TRRUST = '../network-data/TRRUST.csv'
DATA_TFLINK = '../network-data/TFLink.csv'
DATA_STRING = '../network-data/STRING_by_gene.csv'

TF = 'TF'
TARGET = 'Target'


with open(DATA_EXPRESSIONS, 'r') as f_expressions:
    expressions = pd.read_csv(f_expressions)

# with open(DATA_TRRUST, 'r') as f_trrust:
#     trrust = pd.read_csv(f_trrust)

def read_dataframe_subset(path, filter_list):
    with open(path, 'r') as f:
        complete_df = pd.read_csv(f)

    return complete_df[(complete_df[TF].isin(filter_list))
                       & complete_df[TARGET].isin(filter_list)]

def get_TFLink_subset(filter_list):
    return read_dataframe_subset(DATA_TFLINK, filter_list)[['TF', 'Target']]

def get_STRING_subset(filter_list):
    return read_dataframe_subset(DATA_STRING, filter_list)[['TF', 'Target', 'combined_score']]


class WeightMethod:
    DEFAULT_WEIGHT_NO_STRING = 0.1

    def __init__(self, deg, STRING_gene_data):
        self.deg = deg
        self.df = STRING_gene_data

    def GM(self, gene_a, gene_b):
        return sqrt(self.deg[gene_a] * self.deg[gene_b])

    def RMS(self, gene_a, gene_b):
        return sqrt((self.deg[gene_a]**2 + self.deg[gene_b]**2) / 2)

    def STRING(self, gene_a, gene_b):
        condition = ((self.df['gene1'] == gene_a)
                     & (self.df['gene2'] == gene_b))
        result = self.df.loc[condition, 'combined_score']
        if not result.empty:
            # n.b. result.iloc[0] gives you STRING score in thousandths
            return result.iloc[0] / 1000
        else:
            # n.b. arbitrary rating for case where no STRING match is found
            return self.DEFAULT_WEIGHT_NO_STRING

    def product_STRING(self, gene_a, gene_b):
        # Use STRING weight and multiply by RMS
        return self.RMS(gene_a, gene_b) * self.STRING(gene_a, gene_b)


class WeightVectorMethod(WeightMethod):
    DEFAULT_WEIGHT_NO_STRING = 0.1

    def RMS(self, gene, neighbours):
        neighbour_expressions = np.array([self.deg[n] for n in neighbours[TARGET]])
        return np.sqrt((self.deg[gene]**2 + neighbour_expressions**2) / (len(neighbours) + 1))

    def STRING(self, gene, neighbours):
        condition = (self.df[TF] == gene)
        result = self.df.loc[condition, 'combined_score']
        if not result.empty:
            # n.b. result.iloc[0] gives you STRING score in thousandths
            return result / 1000
        else:
            # n.b. arbitrary rating for case where no STRING match is found
            return self.DEFAULT_WEIGHT_NO_STRING

    def product_STRING(self, gene, neighbours):
        return self.RMS(gene, neighbours) * self.STRING(gene, neighbours)
