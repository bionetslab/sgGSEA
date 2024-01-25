import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from itertools import combinations
import os 

# @LIZA: put in the path:
# in your case it will look something like this "C:\Users\annam\sgGSEA\venn_diagrams" (the backslashes on windows need to be from top left to bottom right)
data_path = "/data/bionets/je30bery/hiwi/sgGSEA/venn_diagrams"

# @ LIZA: if you only want to produce one file:
# 1) put Hashtags in front of the for loops (the three lines below that start with "for") and the line after that 
# 2) remove the hashtags of the first three lines (that define the variables) to uncomment them and insert your values (the name variable needs to be the same spelling as here: "T cell migration", "Inflammation", "Cell death", "Proteinopathy", "Neuroimmune interaction"
# if you want to produce all files: leave it like this 

def main():
    # alpha = 0.05
    # condition = "pd"
    # name = "T cell migration"
    # retrieve_venn_diagram_from_csv(condition, name, alpha)
    
    for condition in ["pd", "ibd"]:
        for alpha in [0.1, 0.2, 0.3, 0.4]:
            for name in ["T cell migration", "Inflammation", "Cell death", "Proteinopathy", "Neuroimmune interaction"]:
                retrieve_venn_diagram_from_csv(condition, name, alpha)
    


def retrieve_venn_diagram_from_csv(condition, name, alpha):
    result_csv = pd.read_csv(os.path.join(data_path, condition, f"{condition}_{name}_alpha={alpha}.csv"))
    result_csv.rename({"Unnamed: 0": "Gene"}, inplace=True, axis=1)
    result_csv.set_index("Gene", inplace=True)
    
    # get all possible combinations with 3 query gene sets:
    res_3_combinations = list(combinations(result_csv.columns, 3))

    # @LIZA: here you can change the theme as described here https://seaborn.pydata.org/generated/seaborn.set_theme.html

    # @LIZA: here you can change the colors as described here https://seaborn.pydata.org/generated/seaborn.color_palette.html
    # if you need any specific colors, ChatGPT knows how to put them into a seaborn palette :D

    sns.set_theme()
    pal = sns.color_palette("magma", 4)

    
    if condition == "pd":
        idxs_3 = list(combinations([0,1,2,3], 3))
        
        # @LIZA here you can change the figure size
        # if you want everything to be in a row or a column, replace (2, 2) by (4, 1) or (1, 4)
        fig, axs = plt.subplots(2, 2, figsize=(20,10))

    elif condition == "ibd":
        idxs_3 = list(combinations([0,1,2], 3))
        
        # @LIZA here you can change the figure size
        fig, axs = plt.subplots(1, 1, figsize=(10,10))
    
    for i in range(len(res_3_combinations)):            
        vals_a = result_csv[res_3_combinations[i][0]][result_csv[res_3_combinations[i][0]] < alpha].index
        vals_b = result_csv[res_3_combinations[i][1]][result_csv[res_3_combinations[i][1]] < alpha].index
        vals_c = result_csv[res_3_combinations[i][2]][result_csv[res_3_combinations[i][2]] < alpha].index
        
        string_cutoff = len('FDR corrected P-value ')
        set_labels = [k[string_cutoff:] for k in res_3_combinations[i]]
        if condition == "pd":
            col = int(i/2) 
            row = i % 2
            # @LIZA if you want everything in one column/row you need to replace axs[col, row] by axs[i]
            ax = axs[col, row]
        elif condition == "ibd":
            ax = axs
        venn3([set(vals_a), set(vals_b), set(vals_c)], set_labels=set_labels, ax=ax, set_colors=[pal[idxs_3[i][0]], pal[idxs_3[i][1]], pal[idxs_3[i][2]]], alpha=0.8)

    # Here you can change the title
    fig.suptitle(f"Venn Diagrams: {condition.upper()}: {name}")
    plt.savefig(os.path.join(data_path, condition, f"{condition}_{name}_alpha={alpha}.pdf"))


if __name__ == "__main__":
    main()
