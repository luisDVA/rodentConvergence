# R Code and data for "Prevalence and patterns of convergent ecomorphological evolution in rodents"

Supplementary data and code for reproducing all results and figures

### R 
This folder contains R code for all data analysis and visualization. File names group different stages together (pre-processing, analysis, post-processing, visualization).

### out 
Results and intermediate outputs from l1ou models and their subsequent iterative processing: 

- `processed_regWs`: wheatsheaf index values for all trees and traits analyzed, grouped by regime 
- `processed_edges`: convergent species pairs for each trait for creating the network graphs
- `processed_l1ouModels`: l1ou model objects for each trait run for a block of 100 trees each, in .rds format

### data
Morphological, taxonomic, and phylogenetic data used in the analyses. Raw specimen data is provided for craniodental and mandibular measurements in **cranial_out.csv**, for external measurements in **external_out.csv**, and all specimen details can be found in **specimens_crosswalk.csv**.

----

This dataset ('R Code and data for "Prevalence and patterns of convergent ecomorphological evolution in rodents"') is made available under the Open Data Commons Attribution License: http://opendatacommons.org/licenses/by/1.0.
