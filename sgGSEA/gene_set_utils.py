import gseapy as gp


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
