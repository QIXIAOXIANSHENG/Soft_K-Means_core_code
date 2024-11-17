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
sys.path.append(r"/old_Users/zchen190/soft kmeans/run file/090424")
import sk_function_try as skm
names=locals()

path="/old_Users/zchen190/soft kmeans/run file/090424/ygob_lbt001_all_4_species_with_orth_nogene.csv"
reward=[0,0.001,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.025,0.03,0.035,0.04,0.045,0.5,0.6,0.7,0.8,0.9,1,1.2,1.5,2,2.5,3]
lmax_iter=300 #same as default for k-means
lmax_centerdiff=0.01
seed=2
part_path="/old_Users/zchen190/soft kmeans"
standardize=True
k_max=12  # minimum 3
kmethod=True
for i in reward:
    for j in range(3,k_max+1):

        names["samp"+str(i)+"_"+str(j)]=skm.skmeans(path,i,lmax_iter,lmax_centerdiff,seed,part_path,standardize,k_max,kmethod)
            
        names["samp"+str(i)+"_"+str(j)].GetData()
        #names["samp"+str(i)].InitialCNum()# this function has additional input, don't run with other code
        names["samp"+str(i)+"_"+str(j)].c_num=list(range(1,j+1))
        names["samp"+str(i)+"_"+str(j)].InitialCenter()
        #names["samp"+str(i)+"_"+str(j)].GetAllPartition([1,2,3,4,5,6])

if __name__ == "__main__": 
    num_of_process=35
    list_size=128

    for a in reward:
        for k in range(3,k_max+1):

            names["samp"+str(a)+"_"+str(k)].GetStartTime()
            n=names["samp"+str(a)+"_"+str(k)].GetOrthGroupNum()
            l=list(np.unique(names["samp"+str(a)+"_"+str(k)].sample_data["orth_group"]))
            

            #algorithm part
            for j in range(lmax_iter):
            #step 0,define old_center
                names["samp"+str(a)+"_"+str(k)].renew_center_list()
                
                #step 1, define process
                for i in range(1,num_of_process+1):
                    names["p"+str(i)]=multiprocessing.Process(target=names["samp"+str(a)+"_"+str(k)].run_for_all_orth, args=((l[(i-1)*list_size:i*list_size],))) 
                names["p"+str(36)]=multiprocessing.Process(\
                    target=names["samp"+str(a)+"_"+str(k)].run_for_all_orth, args=((l[4480:4523],))) 
                #  names["p"+str(36)] args=args=((l[num_of_process*list_size:maximum number of orthogroups],))

                #start process
                for i in range(1,num_of_process+1):
                    names["p"+str(i)].start()
                names["p"+str(36)].start()
                #join process
                for i in range(1,num_of_process+1):
                    names["p"+str(i)].join()
                names["p"+str(36)].join()
                #step 2, update new_center
                names["samp"+str(a)+"_"+str(k)].update_center()
                
                #step 3, decide stop iter
                if names["samp"+str(a)+"_"+str(k)].EstStopIter():
                    print("iter time:",names["samp"+str(a)+"_"+str(k)].iter+1)
                    break

                #step 4, class iteration time add 1
                if names["samp"+str(a)+"_"+str(k)].iter<names["samp"+str(a)+"_"+str(k)].max_iter-1:
                    names["samp"+str(a)+"_"+str(k)].iterplus()

            names["samp"+str(a)+"_"+str(k)].runlog["iteration time"]=names["samp"+str(a)+"_"+str(k)].iter

            # get final output
            for i in range(1,num_of_process+1):
                names["p"+str(i)]=multiprocessing.Process(target=names["samp"+str(a)+"_"+str(k)].GetFinalGroup_pre, args=((l[(i-1)*list_size:i*list_size],))) 
            names["p"+str(36)]=multiprocessing.Process(target=names["samp"+str(a)+"_"+str(k)].GetFinalGroup_pre, args=((l[4725:4841],))) 
            
            for i in range(1,num_of_process+1):
                names["p"+str(i)].start()
            names["p"+str(36)].start()
                #join process
            for i in range(1,num_of_process+1):
                names["p"+str(i)].join()
            names["p"+str(36)].join()

            names["samp"+str(a)+"_"+str(k)].GetFinalGroup()

            #get total run time
            names["samp"+str(a)+"_"+str(k)].GetTotalTime()

            #save runlog
            savedir="/old_Users/zchen190/soft kmeans/run file/090424/k="+str(k)
            if not os.path.exists(savedir):
                os.makedirs(savedir)
            names["samp"+str(a)+"_"+str(k)].Save_runlog(savedir)
