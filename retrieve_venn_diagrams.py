from tqdm import tqdm
import sgGSEA as sg
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
import os 
from itertools import chain
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
from itertools import combinations


# provided by Eliza:
#t_cell_migration = ["GO:0006935", "GO:0050900", "GO:0050901","GO:0042110", "GO:0030036", "GO:0038032", "GO:0007155", "GO:0007229",  "GO:0030335", "GO:2000408", "GO:0071347", "GO:0050776", "GO:0008369", "GO:0051493"]
#inflammation = ["GO:0006954", "GO:0050729", "GO:0050728", "GO:0071347", "GO:0050727", "GO:0034341", "GO:0002456", "GO:0050871", "GO:0019221"]
#cell_death = ["GO:0012501", "GO:0006915", "GO:0043065", "GO:0043066", "GO:0070265", "GO:0043068", "GO:0006974", "GO:0010941", "GO:0097193", "GO:0097194", "GO:0006979", "GO:2001237", "GO:0001836", "GO:0010506", "GO:0034976", "GO:0034612"]
#proteinopathy = ["GO:0006457", "GO:0051087", "GO:0006517", "GO:0051788", "GO:0050821", "GO:0033554", "GO:0006914", "GO:0051260", "GO:0051787", "GO:0051259", "GO:0071216", "GO:0031647", "GO:0097033", "GO:0050808", "GO:0097035", "GO:1903900"]
#neuroimmune_interaction = ["GO:0007267", "GO:0150078", "GO:0019221", "GO:0038182", "GO:0010976", "GO:0043524", "GO:0030182", "GO:0043525", "GO:1901215", "GO:0071374"]

protein_aggregation = ["GO:0051087"]
autophagy = ["GO:0006914"]
response_to_misfolded_protein = ["GO:0051788"]
cellular_response_to_protein_misfolding = ["GO:0071216"]
protein_deposition = ["GO:0097033"]
data = "data/adj_pvals_only/"

processes = {"Protein aggregation": protein_aggregation,
             "Autophagy": autophagy,
             "Response to misfolded protein": response_to_misfolded_protein,
             "Cellular response to protein misfolding": cellular_response_to_protein_misfolding,
             "Protein deposition": protein_deposition}

#                   "T cell migration": t_cell_migration,
#                 "Inflammation": inflammation,
#                 "Cell death": cell_death,
#                 "Proteinopathy": proteinopathy,
#                 "Neuroimmune interaction": neuroimmune_interaction}


def plot_venn_diagram(condition, name, alpha, results):
    # draw Venn diagrams of all possible combinations with 3 query gene sets:
    res_3_combinations = list(combinations(results.keys(), 3))
    idxs_3 = list(combinations([0,1,2,3], 3))
    
    sns.set_theme()
    pal = sns.color_palette("magma", 4)
    if condition == "pd":
        fig, axs = plt.subplots(2, 2, figsize=(20,10))
    elif condition == "ibd":
        fig, axs = plt.subplots(1, 1, figsize=(10,10))

    for i in range(len(res_3_combinations)):
        col = int(i/2) 
        row = i % 2

        vals_a = results[res_3_combinations[i][0]][results[res_3_combinations[i][0]]["fdr_corrected_p_value"] < alpha]["symbol"].values
        vals_b = results[res_3_combinations[i][1]][results[res_3_combinations[i][1]]["fdr_corrected_p_value"] < alpha]["symbol"].values
        vals_c = results[res_3_combinations[i][2]][results[res_3_combinations[i][2]]["fdr_corrected_p_value"] < alpha]["symbol"].values
        set_labels = [" vs.\n".join(" ".join(k.split("_")[1:-1]).split("vs")) for k in res_3_combinations[i]]
        if condition == "pd":
            ax = axs[col, row]
        elif condition == "ibd":
            ax = axs
        venn3([set(vals_a), set(vals_b), set(vals_c)], set_labels=set_labels, ax=ax, set_colors=[pal[idxs_3[i][0]], pal[idxs_3[i][1]], pal[idxs_3[i][2]]], alpha=0.8)
    
    fig.suptitle(f"Venn Diagrams: {condition.upper()}: {name}")
    plt.savefig(f"/data/bionets/je30bery/hiwi/sgGSEA/venn_diagrams/{condition}_{name}_alpha={alpha}.pdf")


def create_csv_files(condition, name, alpha, results):
    symbols = np.unique(np.concatenate([np.array(results[k]["symbol"]) for k in results]))
    result_df = pd.DataFrame(index=symbols)
    # write csv:
    set_labels = [" vs. ".join(" ".join(k.split("_")[1:-1]).split("vs")) for k in results.keys()]
    set_labels = ["FDR corrected P-value " + l for l in set_labels]
    
    result_df = pd.DataFrame(index=symbols, columns=set_labels)
    for i, k in enumerate(results.keys()):
        # result_df[set_labels[i]] = [(s in results[k]["symbol"].values) for s in result_df.index]      
        for s in symbols:
            try:
                result_df[set_labels[i]].loc[s] = results[k][results[k]["symbol"] == s]["fdr_corrected_p_value"].values[0]
            except IndexError:
                result_df[set_labels[i]].loc[s] = np.nan
    result_df.to_csv(f"/data/bionets/je30bery/hiwi/sgGSEA/venn_diagrams/{condition}_{name}_alpha={alpha}.csv")
    return result_df


def main():
    num_permutations = 1000
    for condition in ["pd", "ibd"]:
        print("+++++++++++++++++++++++++++++", condition.upper(), "++++++++++++++++++++++++++++++++")
        for alpha in tqdm([0.01, 0.05, 0.1, 0.2, 0.3, 0.4]):
            # get query gene sets for condition:   
            filenames = [f for f in os.listdir(data) if f.startswith(condition)]
            query_gene_sets = [sg.get_query_gene_set(os.path.join(data, f), alpha=alpha) for f in os.listdir(data) if f.startswith(condition)]

            for name in tqdm(processes, leave=False):
                # get genes involved in the pathways of the described process
                target_gene_set = sg.get_genes_for_multiple_pathway_process(processes[name])
                # calculate page rank score and (corrected) pvalues
                results = {}            
                for i, qgs in tqdm(enumerate(query_gene_sets)):
                    result = sg.rank_genes('networks/iid_brain_ppi.txt', target_gene_set=target_gene_set, centrality='pagerank', query_gene_set=qgs, sep=",", num_permutations=num_permutations)
                    #results[os.path.splitext(filenames[i])[0]] = result[result["p_value"] < alpha].sort_values("p_value", ascending=True)
                    results[os.path.splitext(filenames[i])[0]] = result.sort_values("fdr_corrected_p_value", ascending=True)
                
                # plot_venn_diagram(condition, name, alpha, results)
                create_csv_files(condition, name, alpha, results)
                
                
if __name__ == "__main__":
    main()