# RNAseq-undergrad-thesis
RNA seq project. Use Soft K-means to cluster RNA seq data of four yeast species based on orthology and expression pattern.

Define Soft K-means class in skm_class.py.

Input data format should have the following properties:
  1. Data should be .csv file.
  2. First two column names should be "orth_group" and "num". The first column is the orthogroup information (integer) for each gene, and the second colunm is the order of the genes in each orthogroup (integer), respectively. Data should be ordered by orthogroup (small to large) and num (small to large, continuous from 1 to ni, where i is the corresponding orthogroup and ni is the size of that orthogroup).
  3. Starting from the third column is the expression value of each gene. Can be named anything.
