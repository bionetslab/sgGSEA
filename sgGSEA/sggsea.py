import graph_tool as gt
import graph_tool.centrality as gtc
import graph_tool.generation as gtg
import pandas as pd


def rank_genes(path_to_network: str, target_gene_set: list[str], centrality: str, sep: str = ' ',
               num_permutations: int = 100, damping: float=0.85, query_gene_set=None) \
        -> pd.DataFrame:
    """Ranks the genes in the query gene set via empirical P-values based on their network centralities w.r.t. the
    target gene set.

    :param path_to_network: path to gene network given as edge list with gene symbols as node IDs; first row is skipped
    :type path_to_network: string
    :param target_gene_set: set of target genes (pathway) given as gene symbols
    :type target_gene_set: list of strings
    :param centrality: centrality measure to be used; select from { 'betweenness', 'pagerank' }
    :type centrality: string
    :param sep: character used as separator in the input network; default: ' '
    :type sep: string
    :param num_permutations: number of times input network is randomized; determines p-value resolution; default: 100
    :type num_permutations: int
    :param damping: damping factor used for pageranke centrality; default: 0.85
    :type damping: float
    :param query_gene_set: list of genes that should be ranked; if None, all genes in the network are ranked
    :type query_gene_set: list of strings or None
    :return: data frame sorted by p-value containing the results for all genes in the query set
    :rtype: pandas.DataFrame
    """
    network, symbol_to_id = _load_network(path_to_network, sep=sep)
    target_gene_set = _prune_gene_set(target_gene_set, symbol_to_id)
    centralities = _compute_centralities(network, target_gene_set, centrality, damping)
    p_values = network.new_vp('double', val=1)
    for _ in range(num_permutations):
        gtg.random_rewire(network)
        centralities_permuted_network = _compute_centralities(network, target_gene_set, centrality, damping)
        p_values.a += (centralities_permuted_network.a >= centralities.a)
    p_values.a /= (num_permutations + 1)
    if query_gene_set is None:
        query_gene_set = list(range(network.num_vertices()))
    else:
        query_gene_set = _prune_gene_set(query_gene_set, symbol_to_id)
    results = pd.DataFrame(data={
        'symbol': [network.vp['name'][v] for v in query_gene_set],
        centrality: [centralities[v] for v in query_gene_set],
        'p_value': [p_values[v] for v in query_gene_set],
        'degree': [network.degree_property_map('total')[v] for v in query_gene_set],
        'in_target_gene_set': [v in target_gene_set for v in query_gene_set]
    })
    results.sort_values(by='p_value', inplace=True)
    results.reset_index(drop=True, inplace=True)
    return results


def _compute_centralities(network: gt.Graph, target_gene_set: list[str],
                          centrality: str, damping: float=0.85) -> gt.VertexPropertyMap:
    if centrality == 'betweenness':
        centralities, _ = gtc.betweenness(network, pivots=target_gene_set)
    elif centrality == 'pagerank':
        p = network.new_vp('double')
        for gene in target_gene_set:
            p[gene] = 1.0
        centralities = gtc.pagerank(network, damping=damping, pers=p)
    return centralities


def _load_network(path_to_network: str, sep: str=' ') -> gt.Graph:
    network = gt.load_graph_from_csv(path_to_network, skip_first=True, csv_options={'delimiter': sep})
    symbol_to_id = {network.vp['name'][v]: int(v) for v in network.vertices()}
    return network, symbol_to_id


def _prune_gene_set(gene_set: list[str], symbol_to_id: dict[str, int]) -> list[str]:
    pruned_gene_set = list(set(gene_set).intersection(set(symbol_to_id.keys())))
    return [symbol_to_id[gene] for gene in pruned_gene_set]
