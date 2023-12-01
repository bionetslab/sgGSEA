import gseapy as gp
import numpy as np

def show_library_names(organism: str = 'human') -> list[str]:
    """Shows names of all queryable gene set libraries.

    :param organism: select one from { 'human', 'mouse', 'yeast', 'fly', 'fish', 'worm' }
    :type organism: string
    :return: list of queryable libraries
    :rtype: list of strings
    """
    return gp.get_library_name(organism=organism)


def get_library(library_name: str, organism: str = 'human') -> dict[str, list[str]]:
    """Downloads gene set library.

    :param library_name: name of gene set library
    :type library_name: string
    :param organism: select one from { 'human', 'mouse', 'yeast', 'fly', 'fish', 'worm' }
    :type organism: string
    :return: gene sets contained in the library
    :rtype: dictionary
    """
    return gp.parser.download_library(name=library_name, organism=organism)


def show_gene_set_names(library: dict[str, list[str]]) -> list[str]:
    return list(library.keys())


def get_gene_set(library: dict[str, list[str]], gene_set_name: str) -> list[str]:
    return library[gene_set_name]

def get_genes_for_multiple_pathway_process(pathways):
    pathways = {p: None for p in pathways}
    for version in ["2013", "2015", "2017", "2017b", "2018", "2021", "2023"][::-1]:
        go_bp = get_library(f'GO_Biological_Process_{version}')
        gene_sets = show_gene_set_names(go_bp)
        for p in pathways:
            if not pathways[p]:
                for g in gene_sets:
                    if p in g:
                        pathways[p] = get_gene_set(go_bp, g)
    
    target_gene_set = list()
    for k in pathways:
        try:
            target_gene_set += pathways[k]
        except:
            print(k, "not found")
    target_gene_set = np.unique(target_gene_set)
    return target_gene_set
        