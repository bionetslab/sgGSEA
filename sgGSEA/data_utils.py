import pandas as pd 
import numpy as np

def get_query_gene_set(path: str, alpha: float = 0.05):
    query_gene_set = pd.read_csv(path)
    query_gene_set = query_gene_set[query_gene_set["padj"] < alpha]
    return query_gene_set["gene_name"]


def df_from_xls(path:str, has_colnames: bool = True):
    with open(path, 'r') as f:
        lines = f.readlines()
    cols = None
    if has_colnames:
        cols = lines[0].split("\t")
        cols[-1] = cols[-1].split("\n")[0]
        lines = lines[1:]
    
    split_lines = list()    
    for i, l in enumerate(lines):
        spl = l.split("\t")
        spl[-1] = spl[-1].split("\n")[0]
        split_lines.append(spl)
        
    df = pd.DataFrame(split_lines, columns=cols)
    return df


def fpkm_expression_matrix(path):
    df = df_from_xls(path)
    df = df.set_index("gene_id")
    fpkms = [c for c in df.columns if "fpkm" in c]
    expression_df = df[fpkms]
    return expression_df