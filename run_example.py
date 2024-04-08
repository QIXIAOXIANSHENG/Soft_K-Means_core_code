import numpy as np
import os
import copy
import pandas as pd
from itertools import combinations
from sklearn.preprocessing import MinMaxScaler
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import sys
import datetime
import pickle
import multiprocessing

# this path is the file that contains the .py file defining soft k-means class
sys.path.append(r"/old_Users/zchen190/soft kmeans/run file/040724")
# import the .py file defining soft k-means class
import sk_function_try as skm

names=locals()

# input RNA-seq data path
path="/old_Users/zchen190/soft kmeans/all_ygob_data/all_4_species_nogene.csv"

# output path of the algorithm
op_path="/old_Users/zchen190/soft kmeans/run file/040724"

# define different rewards
reward=[0.001,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.025,0.03,0.035,0.04,0.045]

# maximum iterations
lmax_iter=300 #same as default for k-means

# stop signal, maximum distance between centers to stop iteration
lmax_centerdiff=0.01

# seed for k-means initializaion
seed=2

# path to store partition information
part_path="/old_Users/zchen190/soft kmeans"

# standardization (0-1) of input RNA-seq data
standardize=True

# maximum number of clusters
k_max=10

# use k-means to initialize
kmethod=True

# define classes with different reward w
for i in reward:
    names["samp"+str(i)]=skm.skmeans(path,i,lmax_iter,lmax_centerdiff,seed,part_path,standardize,k_max,kmethod)
        
    names["samp"+str(i)].GetData()
    # this would output elbow plot for k-means, and user need to input the number of clusters
    names["samp"+str(i)].InitialCNum()# this function has additional input, don't run with other code
    # or define number of clusters manually. There are 8 clusters in this example
    names["samp"+str(i)].c_num=[1,2,3,4,5,6,7,8]
    names["samp"+str(i)].InitialCenter()
    # every possible partition of different orthogroup sizes. The example has maximum orthogroup size 6
    names["samp"+str(i)].GetAllPartition([1,2,3,4,5,6])

# run for different classes

if __name__ == "__main__": 
	# should check for number of orthogroups before running code to define appropriate process number and data used in each process
    num_of_process=35
    list_size=141 # each process in the first 35 processes calculates 141 orthogroups

    for a in reward:

        names["samp"+str(a)].GetStartTime()
        n=names["samp"+str(a)].GetOrthGroupNum()
        l=list(np.unique(names["samp"+str(a)].sample_data["orth_group"]))
        
        #num_of_process=10
        #list_size=100

        

        #algorithm part
        for j in range(lmax_iter):
        #step 0,define old_center
            names["samp"+str(a)].renew_center_list()
            
            #step 1,
            for i in range(1,num_of_process+1):
                names["p"+str(i)]=multiprocessing.Process(target=names["samp"+str(a)].run_for_all_orth, args=((l[(i-1)*list_size:i*list_size],))) 
		    # define another process because the number of orthogroups could not be divided by 36 (number of nodes used)
            names["p"+str(num_of_process+1)]=multiprocessing.Process(target=names["samp"+str(a)].run_for_all_orth, args=((l[4935:5011],))) 

            #start process
            for i in range(1,num_of_process+1):
                names["p"+str(i)].start()
            names["p"+str(num_of_process+1)].start()
            #join process
            for i in range(1,num_of_process+1):
                names["p"+str(i)].join()
            names["p"+str(num_of_process+1)].join()
            #step 2, update new_center
            names["samp"+str(a)].update_center()
            
            #step 3, decide stop iter
            if names["samp"+str(a)].EstStopIter():
                print("iter time:",names["samp"+str(a)].iter+1)
                break

            #step 4, class iteration time add 1
            if names["samp"+str(a)].iter<names["samp"+str(a)].max_iter-1:
                names["samp"+str(a)].iterplus()

        names["samp"+str(a)].runlog["iteration time"]=names["samp"+str(a)].iter

        # get final output
        for i in range(1,num_of_process+1):
            names["p"+str(i)]=multiprocessing.Process(target=names["samp"+str(a)].GetFinalGroup_pre, args=((l[(i-1)*list_size:i*list_size],))) 
        names["p"+str(num_of_process+1)]=multiprocessing.Process(target=names["samp"+str(a)].GetFinalGroup_pre, args=((l[4935:5011],))) 
        
        for i in range(1,num_of_process+1):
            names["p"+str(i)].start()
        names["p"+str(num_of_process+1)].start()
            #join process
        for i in range(1,num_of_process+1):
            names["p"+str(i)].join()
        names["p"+str(num_of_process+1)].join()

        names["samp"+str(a)].GetFinalGroup()

        #get total run time
        names["samp"+str(a)].GetTotalTime()

        #save runlog, parameter is the path to save the output log
        names["samp"+str(a)].Save_runlog(op_path)