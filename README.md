Floral symmetry, visitation, and sharing
========================================

This repository contains data and scripts supporting 

Yoder JB, G Gomez, and CJ Carlson. Zygomorphic flowers have fewer visitors. *bioRxiv*, doi: [10.1101/743872](https://www.doi.org/10.1101/743872)

Floral visitation data are derived from records archived in the [Web of Life](http://www.web-of-life.es) and [Interaction Web DataBase](https://www.nceas.ucsb.edu/interactionweb) repositories; source datasets and the repositories where each were obtained are given in [`data/references_all.csv`](data/references_all.csv).

Updated as of 3 July 2020.


Contents
--------

### data

Files describing original source material

- `references_all.csv` --- table of source datasets, metadata, and original references

### output

Files containing data derived from compiling source data

- `plant_degree-sharing.csv` --- floral symmetry, visitor count, and visitor sharing for all plant species 
- `whole-matrix_stats.csv`--- network structure statistics for all plant-visitation networks in the dataset
- `sub-matrix_stats.csv` --- network structure statistics for all plant-visitation networks, divided into sub-networks based on plant floral symmetry
- `degree_per_plant_phy.csv` --- floral symmetry, visitor count, and visitor sharing for all plant species, with principal components of the phylogenetic distances among them, based on the supertree constructed by [Smith and Brown (2018)](https://doi.org/10.1002/ajb2.1019)
- `figures` --- sub-folder containing all figures

### scripts

- `all-data_processing.R` --- Generating working files for further analysis
	- requires: network data, `plant_degree-sharing.csv`, `whole-matrix_stats.csv`, `sub-matrix_stats.csv`, `references_all.csv`
	- generates: `full-network_stats-meta.txt`, `sub-network_stats-meta.txt`, `plant_stats-meta.txt`
- `plant_degree_analysis.R` --- species-level analysis, model-fitting
	- requires: `degree_per_plant_phy.csv`, `references_all.csv`
	- generates: model-fitting and comparison results
- `figures.R` --- sub-network analyses, figure generation
	- requires: `full-network_stats-meta.txt`, `sub-network_stats-meta.txt`, `plant_stats-meta.txt`
	- generates: `Fig01_map-symm-npoll-Rshare.pdf`, `Fig02_pairs-cors.pdf`
