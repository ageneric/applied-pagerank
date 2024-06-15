import rbo
import scipy


def compute_kl_divergence(p, q):
    """p, q: PageRank vectors."""
    return sum(scipy.special.rel_entr(p, q))

def compute_rbo(order, other_order, k=None, p=1):
    """order, other_order: lists."""
    return rbo.RankingSimilarity(order, other_order).rbo(k, p)

def compute_rbo_less_equal_prs(order, other_order, k=None, p=1):
    """order, other_order: lists of tuples."""

    # Find modal pagerank of order and other_order
    # note: assumes pageranks are not multi-modal
    mode1 = max([x[1] for x in order], key=[x[1] for x in order].count)
    mode2 = max([x[1] for x in other_order], key=[x[1] for x in other_order].count)

    # Remove items with modal pagerank from order and other_order
    list1 = [x[0] for x in order if x[1] != mode1]
    list2 = [x[0] for x in other_order if x[1] != mode2]
    return rbo.RankingSimilarity(list1, list2).rbo(k, p)
