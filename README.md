# Soft K-Means
RNA seq project. Use Soft K-means to cluster RNA seq data of four yeast species based on orthology and expression pattern.
The Soft K-means is defined in the SKM_define branch. Last modified on 4/8/24.

# Differentially expressed orthologous genes detection pipeline
1) Gene Selection: Identifying genes with valid expression patterns under stress, using time series analysis; 
2) Co-Clustering: Clustering genes based on expression patterns and prior orthology information, guided by a set of algorithm parameters; 
3) Grid Search: Conducting a search to find the optimal clustering within the cluster space; 
4) Differentially Expressed Orthologous Gene Reporting: Labeling and reporting orthologous genes that fall into clusters differing from the majority. 

# Use of code
Use the codes in the following order:
  gene_selection.R
  merge_selected_genes.R
  ygob_generate.R
  find_optimal_K.py
  run_argon.sh
  generate_cluster_plotcsv.py
  plot_heatmap.R
  statistical_scores.py
  topGO_enrichment_analysis.R
  biological_scores.R
  combind_score.R

# Soft K-Means input data format
Input data format should have the following properties:

1. Data should be .csv file.
2. First two column names should be "orth_group" and "num". The first column is the orthogroup information (integer) for each gene, and the second colunm is the order of the genes in each orthogroup (integer), respectively. Data should be ordered by orthogroup (small to large) and num (small to large, continuous from 1 to ni, where i is the corresponding orthogroup and ni is the size of that orthogroup).
3. Starting from the third column is the expression value of each gene. Can be named anything.

# Soft K-Means output data format
The output pickle file contains a Python dictionary with the used data, centers in each iteration, the final cluster assign, and the run time.
