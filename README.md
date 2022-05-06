# Spectral Clustering Analysis of GDSC Data

A brief video summary of this project can be found [here](https://youtu.be/J3tc1sJzkm0) on YouTube. 

## /assets

Contains relevant csv files

### Raw data

+ **GDSC2_fitted_dose_response_25Feb20.csv** - contains the raw IC<sub>50</sub> data for 135,242 drug / cell line combinations; data can be downloaded [here](ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_25Feb20.xlsx)
+ **cell_list.csv** - contains GDSC-generated information about the included cell lines, including tissue and TCGA classification; data can be downloaded [here](ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx)
+ **drug_data.csv** - contains GDSC-generated information about the included drugs, including drug pathways and targets; data can be found [here](https://www.cancerrxgene.org/downloads/drug_data)

## /src

* **main.R** - contains
* **main.ipynb** - contains 
* **graph.py**
* **kmeans.py**
    * `find_kmeans`: find an optimal number clusters via the elbow method and fit *k*-means with this many clusters
    * `plot_kmeans`: plot the SSE vs. clusters and elbow point for cell and drug data
* **utils.py**
    * `import_data`: loads GDSC data and pre-process into a wide matrix
    * `process_data`: produces mean-centered data and masks