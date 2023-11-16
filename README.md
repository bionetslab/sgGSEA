# sgGSEA
Network-based single-gene gene set enrichment analysis

# Dependencies

sgGSEA requires gseapy and graph-tool, as well as jupyter and seaborn (to run tutorial.py). 

The following workflow should work to install these dependencies:

```
conda create --name sggsea -c conda-forge graph-tool
conda activate sggsea
conda install pip
pip install gseapy
pip install jupyter
pip install seaborn
python -m ipykernel install --user --name=sggsea
```


