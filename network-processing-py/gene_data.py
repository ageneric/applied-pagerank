def get_gene_differential_expressions(expressions):
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
