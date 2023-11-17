import pandas as pd 

def get_query_gene_set(path: str, alpha: float = 0.05):
    query_gene_set = pd.read_csv(path)
    query_gene_set = query_gene_set[query_gene_set["padj"] < alpha]
    return query_gene_set["gene_name"]