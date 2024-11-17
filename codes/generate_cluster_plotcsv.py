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
import pickle

names=locals()

reward=[0.6,0.7,0.8,0.9,1,1.5]
k_list=[3,4,5,7,9]

for w in reward:
    for k in k_list:
        filename=r"D:\Grad study\lab\data\soft kmeans\gene selection\ygob_libnorm_lbt005_cluster_result\pickle file\k="+str(k)+"\w="+str(w)+".pkl"
        with open(filename,"rb") as f_obj:
           a=pickle.load(f_obj)
        names["r_"+str(k)+"_"+str(w)]=a["final group"]
        names["r_"+str(k)+"_"+str(w)]=names["r_"+str(k)+"_"+str(w)].sort_values(by=["orth_group","num"])
        names["p_"+str(k)+"_"+str(w)]=pd.merge(left=names["r_"+str(k)+"_"+str(w)],right=a["data"],how='left',\
                                               left_on=['orth_group', 'num'],right_on=['orth_group', 'num'])


for w in reward:
    for k in k_list:
        filename_c=r"D:\Grad study\lab\data\soft kmeans\gene selection\ygob_libnorm_lbt005_cluster_result\cluster file\k="+str(k)
        filename_p=r"D:\Grad study\lab\data\soft kmeans\gene selection\ygob_libnorm_lbt005_cluster_result\plot_heatmap_csv file\k="+str(k)
        if not os.path.exists(filename_c):
            os.mkdir(filename_c)
        if not os.path.exists(filename_p):
            os.mkdir(filename_p)
        names["r_"+str(k)+"_"+str(w)].to_csv(filename_c+"\w="+str(w)+".csv",index=False)
        names["p_"+str(k)+"_"+str(w)].to_csv(filename_p+"\w="+str(w)+".csv",index=False)
