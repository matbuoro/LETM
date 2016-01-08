Genetic architecture
====================================


This folder contains the datasets (.RData or .csv)  and the R scripts to run the analysis for the 'half-sib common environment' (named 'ISLANDS') and a 'family level split environment' (named 'EXPERIMENTS') experiments.

Note that we use the population codes :
BR: North Sea islands of Brownsman,
WWO: West Wideopen
EWO: East Wideopen

Overview
--------

Earwigs
	|_ code/           # any programmatic code
	|	
	|_ data            # raw and primary data, are not changed once created 
	|	
	|_ docs/            # documentation for the study
	|  +- article/       # manuscript(s), whether generated or not
	|_ results         # all output from workflows and analyses
	|
	|- README          # the top level description of content


The *code* folder contains the scripts (R, SQL,...) used to clean the dataset and run the analysis. 

The *data* folder contains the datasets used to run the analysis and some scripts to clean up the raw dataset.

All outputs (figures,...) are kept in the *results* folder. 

The *docs* folder contains documents such as presentation, relevant bibliography, ...


