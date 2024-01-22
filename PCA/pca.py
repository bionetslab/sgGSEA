import pandas as pd
from sklearn.decomposition import PCA
import sgGSEA as sg
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import json


# @ELIZA: choose problem: "ibd" or "pd"
problem = "ibd"

# @ELIZA: put in the path:
# in your case it will look something like this f"C:\Users\annam\sgGSEA\data\original\{problem}"
# the "f" needs to be in front of the address, so that the value of the variable "problem" gets inserted
base_path = f"/data/bionets/je30bery/hiwi/sgGSEA/data/original/{problem}"


# load fpkm values
expression_matrix = list()
for f in os.listdir(base_path):
    expression_df = sg.fpkm_expression_matrix(os.path.join(base_path, f))
    expression_matrix.append(expression_df)
expression_matrix = pd.concat(expression_matrix, axis=1)
# drop duplicate samples
expression_matrix = expression_matrix.T.drop_duplicates().T
# drop genes with NAN values in expression -> genes that were not analyzed for all samples
expression_matrix = expression_matrix.dropna()
# perform f(x)=log(x+1), +1 to avoid -np.inf values for zeroes 
expression_matrix_log = np.log1p(expression_matrix.astype(float))
# use PCA
pca = PCA(2)
transformed = pca.fit_transform(expression_matrix_log.T)


# load label mappings
# @ELIZA: if you want to change the names in the legend, you need to open the file "mappings.json" and change them there
with open('./PCA/mappings.json', 'r') as json_file:
    label_mapping  = json.load(json_file)[problem]
    
    
# map labels
expression_matrix.columns = [label_mapping["_".join(c.split("_")[:-1])] for c in expression_matrix.columns]
# plot and color by label, save plot

# @ELIZA: here you can change the theme as described here https://seaborn.pydata.org/generated/seaborn.set_theme.html
sns.set_theme("notebook")

# @ELIZA: here you can change the colors as described here https://seaborn.pydata.org/generated/seaborn.color_palette.html
# if you need any specific colors, ChatGPT knows how to put them into a seaborn palette :D

pal = sns.color_palette("magma", len(np.unique(expression_matrix.columns)))

for i, condition in enumerate(np.unique(expression_matrix.columns)):
    condition_samples = transformed[np.where(expression_matrix.columns == condition)]
    plt.scatter(condition_samples[:,0], condition_samples[:,1], label=condition, color=pal[i])
    
# @ELIZA: here you can change things like fontsize or title or where the file will be saved
plt.legend(fontsize=8)
plt.title(f"PCA on log-transformed expression data of the {problem.upper()} dataset")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(f"pca_{problem}.pdf")

