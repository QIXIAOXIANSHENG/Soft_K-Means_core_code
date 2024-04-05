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

names=locals()

class skmeans():
    def __init__(self,path,reward,lmax_iter,lmax_centerdiff,seed,part_path,standardize=True,k_max=np.nan,kmethod=False):
        self.path=path
        self.w=reward
        self.max_iter=lmax_iter
        self.max_centerdiff=lmax_centerdiff
        self.kmeans_method=kmethod
        self.seed=seed
        self.runlog=multiprocessing.Manager().dict({})
        self.standardize=standardize
        self.part_path=part_path
        self.k_max=k_max
        self.sample_data=pd.DataFrame()
        self.sample_data_gene=[]
        self.c_num=[]
        self.new_center=pd.DataFrame()
        self.old_center=pd.DataFrame()
        self.names=locals()
        self.iter=0
        self.finalgroup=multiprocessing.Manager().list([])
        

        

    def GetData(self):
        df=pd.read_csv(self.path)
        if self.standardize:
            values = df.iloc[:,2:].T.values 
            values = values.astype('float64') 


            tool = MinMaxScaler(feature_range=(0, 1))
            data = tool.fit_transform(values)
    
            self.sample_data=pd.DataFrame(data).T
        
            self.sample_data.insert(0, "orth_group", df["orth_group"])
            self.sample_data.insert(1, "num", df["num"])
            self.sample_data.columns=df.columns
                    
        else:
            self.sample_data=df
        self.sample_data_gene=self.sample_data.columns[2:]

        self.runlog["data"]=self.sample_data


    def GetOrthGroupNum(self):
        return len(list(self.sample_data["orth_group"].value_counts().index))
    
    def InitialCNum(self):
        if self.kmeans_method:
            n=self.GetOrthGroupNum()
            distortions = []
            inertias = []
            mapping1 = {}
            mapping2 = {}

            if np.isnan(self.k_max):
                self.k_max=int(n*1.2)+1
            K = range(1, self.k_max+1)
            X=np.array(self.sample_data[self.sample_data_gene])
            for k in K:
                kmeanModel = KMeans(n_clusters=k).fit(X)
                kmeanModel.fit(X)
                distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_,'euclidean'), axis=1)) / X.shape[0])
                inertias.append(kmeanModel.inertia_)
                mapping1[k] = sum(np.min(cdist(X, kmeanModel.cluster_centers_,'euclidean'), axis=1)) / X.shape[0]
                mapping2[k] = kmeanModel.inertia_
            plt.plot(K, inertias, 'bx-')
            plt.xlabel('Values of K')
            plt.ylabel('Inertia')
            plt.title('The Elbow Method using Inertia')
            plt.show()
            
            k=int(input("choose k:"))
            self.c_num=list(range(1,k+1))
        
        else:
            self.c_num=list(self.sample_data["orth_group"].value_counts().index)
            self.c_num.sort()

        self.runlog["used values"]={"path":self.path,"reward":self.w,"max iter":self.max_iter,
                                    "max center difference":self.max_centerdiff,"seed":self.seed,
                           "standardize":self.standardize,"start center amount":len(self.c_num),
                           "if use k-means":self.kmeans_method}

    def InitialCenter(self):
        np.random.seed(self.seed)
        index=["center_num"]+list(self.sample_data_gene)
        insert=self.c_num.copy()
        
        if self.kmeans_method:
            X=np.array(self.sample_data[self.sample_data_gene])
            kmeanModel = KMeans(n_clusters=len(self.c_num)).fit(X)
            a=kmeanModel.cluster_centers_
        else:
            a=self.sample_data.groupby("orth_group").mean().values[:,1:]+\
                np.random.normal(0,0.05,len(self.sample_data_gene)*len(self.c_num)).reshape(len(self.c_num),len(self.sample_data_gene))

        data=np.insert(a,0,insert,axis=1)
        df=pd.DataFrame(data,columns=index)
        self.new_center=df
        self.runlog["start center"]=df

    def GetOrthGroupLength(self,orth_group):
        return self.sample_data["orth_group"].value_counts()[orth_group]
    
    def GenerateDataSaver(self,orth_group):
        index=["compute","diff_origin"]+["d_"+str(c) for c in self.c_num]+["group","distance"]
        end_list = []
        list1=list(range(1,self.GetOrthGroupLength(orth_group)+1))
        for i in range(1,len(list1)+1):
            end_list.extend(combinations(list1, i))
        colnames=[i for i in end_list]
        return pd.DataFrame(columns=colnames,index=index)
    
    def GetAllPartition(self,list1):
        filename=self.part_path+"/all partition/"+str(1)+".pkl"
        with open(filename, 'wb') as f_obj:
            pickle.dump([[tuple([list1[0]])]], f_obj) 
        
        if len(list1)==1:
            return
        
        for r in range(2,len(list1)+1):
            prev_filename=self.part_path+"/all partition/"+str(r-1)+".pkl"
            with open(prev_filename,"rb") as f_obj:
                prev_part = pickle.load(f_obj)
            new_el=list1[r-1]
            new_part=[]
            for i in prev_part:
                new_part.append(i+[tuple([new_el])])
                for j in i:
                    new=tuple(list(j)+[new_el])
                    el=[k for k in i if k!=j]
                    el.append(new)
                    new_part.append(el)
                    
            new_filename=self.part_path+"/all partition/"+str(r)+".pkl"
            with open(new_filename, 'wb') as f_obj:
                pickle.dump(new_part, f_obj)

    def GetMinimalGroup(self,orth_group):#with penalty for partition num
        data_saver=self.GenerateDataSaver(orth_group)
        def CalculateDistance(length):
            if length==1:
                for i in self.c_num:
                    names["dist_with"+str(i)]=self.sample_data[self.sample_data["orth_group"]==orth_group][self.sample_data_gene]-np.array(self.old_center[self.old_center["center_num"]==i][self.sample_data_gene])
                    names["dist_with"+str(i)]=names["dist_with"+str(i)].apply(func=np.linalg.norm,axis=1)
                df_dist_orth=pd.concat([names["dist_with"+str(i)] for i in self.c_num],axis=1)
                df_dist_orth.columns=["d_"+str(i) for i in self.c_num]
                df_dist_orth.index=[tuple([i]) for i in list(range(1,self.GetOrthGroupLength(orth_group)+1))]
                df_dist_orth=df_dist_orth.T
                data_saver.loc[data_saver.index.str.contains('d_'),tuple((df_dist_orth.columns))]=df_dist_orth
                data_saver.loc[data_saver.index.str.contains('compute'),tuple((df_dist_orth.columns))]=True
                data_saver.loc[data_saver.index.str.contains('diff_origin'),tuple((df_dist_orth.columns))]=1
                b=data_saver.loc[data_saver.index.str.contains('d_'),tuple((df_dist_orth.columns))].apply(min,axis=0)
                
                x=data_saver.loc[data_saver.index.str.contains('d_'),:]
                y=list([x.index[x[i]==b[i]][0] for i in b.index])
                data_saver.loc[data_saver.index.str.contains('group'),tuple((df_dist_orth.columns))]=y
                #data_saver.loc[data_saver.index.str.contains('group'),tuple((df_dist_orth.columns))]=[list(data_saver.index[data_saver[i]==b[i]])[0] for i in b.index]
                data_saver.loc[data_saver.index.str.contains('distance'),tuple((df_dist_orth.columns))]=pd.DataFrame(b,columns=["distance"]).T
            else:
                col_len_is_length=[i for i in data_saver.columns if len(i)==length]
                partition_col=[[col_len_is_length[i][0:(len(col_len_is_length[i])-1)],tuple([col_len_is_length[i][len(col_len_is_length[i])-1]])] for i in range(len(col_len_is_length))]
                for i in range(len(partition_col)):
                    data_origin=[tuple([j]) for j in col_len_is_length[i]]
                    
                    if data_saver.loc[data_saver.index.str.contains('compute')][partition_col[i][0]][0]==True:
                        origin_group_num=data_saver[data_origin].loc["group"].value_counts()
                        #origin_group=data_saver[data_origin].loc["group"]
                        c=data_saver.loc[data_saver.index.str.contains('d_')][partition_col[i][0]]+data_saver.loc[data_saver.index.str.contains('d_')][partition_col[i][1]]
                        data_saver.loc[data_saver.index.str.contains('d_'),[col_len_is_length[i]]]=pd.DataFrame(c,columns=[col_len_is_length[i]])
                        c_copy=c.copy()
                        if len(origin_group_num)>1:
                            data_saver.loc[data_saver.index.str.contains('diff_origin'),[col_len_is_length[i]]]=len(origin_group_num)
                            for k in list(origin_group_num.index):
                                c_copy[k]=c_copy[k]-self.w*(sum(origin_group_num)-origin_group_num[k])
                        else:
                            data_saver.loc[data_saver.index.str.contains('diff_origin'),[col_len_is_length[i]]]=1
                            org_group=list(origin_group_num.index)[0]
                            c_copy[c.index!=org_group]=c[c.index!=org_group]-self.w
                            
                            
                        #group,distance=c_copy.index[c_copy==min(c_copy)],c[c_copy==min(c_copy)][0]
                        group,distance=c_copy.index[c_copy==min(c_copy)],c_copy[c_copy==min(c_copy)][0]
                        data_saver.loc[data_saver.index.str.contains('group'),[col_len_is_length[i]]]=group
                        data_saver.loc[data_saver.index.str.contains('dist'),[col_len_is_length[i]]]=distance
                        data_saver.loc[data_saver.index.str.contains('compute'),[col_len_is_length[i]]]=True
                    else:
                        CalculateDistance(length-1)
                        CalculateDistance(length)
            
        
        length=self.GetOrthGroupLength(orth_group)
        CalculateDistance(length)
        group_assign=[int(i.split("_")[1]) for i in [i for _,i in data_saver.loc["group"].items()]]
        data_saver.loc[data_saver.index.str.contains('group'),data_saver.columns]=pd.DataFrame(group_assign,columns=["group"],index=data_saver.columns).T
        
        filename=self.part_path+"/all partition/"+str(length)+".pkl"
        with open(filename,"rb") as f_obj:
            possible_partition=pickle.load(f_obj)
        
        b=[sum(np.array(data_saver.loc[["distance"],i])[0]) for i in possible_partition]
        c=[b[i]+1e-12*len(possible_partition[i]) for i in range(len(b))]
        mini_index=c.index(min(c))
        final_group_assign=data_saver.loc[["group"],possible_partition[mini_index]]
        return (data_saver,final_group_assign)
    
    def renew_center_list(self):
        self.old_center=self.new_center
        for i in self.c_num:
            self.names["center_"+str(i)]=multiprocessing.Manager().list([])

    def run_for_all_orth(self,list_of_orth):#parallelize!!!!!!!
        for i in list_of_orth:# this list_of_orth contains the orthogroup list itself, not its index anymore.
            a_val=self.GetMinimalGroup(i)
            orth_assign=a_val[1]
            #self.runlog["iter "+str(self.iter)+" orth "+ str(i)+" data"]=a_val
            centers=list(np.array(orth_assign)[0])
            cols=list(orth_assign.columns)
            for j in range(len(centers)):
                for k in list(cols[j]):
                    m=np.array(self.sample_data[self.sample_data["orth_group"]==i][self.sample_data["num"]==k][self.sample_data_gene])[0]
                    self.names["center_"+str(centers[j])].append(m)
                #m=[np.array(self.sample_data[self.sample_data["orth_group"]==i][self.sample_data["num"]==k][self.sample_data_gene])[0] for k in list(cols[j])]
                



    def update_center(self):
        c_num_copy=self.c_num.copy()
        for i in c_num_copy:
            try:
                self.names["center_"+str(i)]=sum(self.names["center_"+str(i)])/len(self.names["center_"+str(i)])
            except:
                self.c_num.remove(i)
                print("loss one center.")

        index=["center_num"]+list(self.sample_data_gene)
        insert=self.c_num.copy()
        a=[self.names["center_"+str(i)][j] for i in self.c_num for j in range(len(self.names["center_"+str(i)]))]
        a=np.array(a).reshape(len(self.c_num),len(self.sample_data_gene))
        data=np.insert(a,0,insert,axis=1)
        self.new_center=pd.DataFrame(data,columns=index)
        self.runlog["iter "+str(self.iter+1)+" center"]=self.new_center
            

    def EstStopIter(self):
        if self.old_center.shape==self.new_center.shape:
            lose_center=False
            diff=np.linalg.norm(np.array(self.old_center[self.sample_data_gene])-np.array(self.new_center[self.sample_data_gene]))
            if diff<=self.max_centerdiff:
                self.runlog["iter "+str(self.iter+1)+" lose center"]=lose_center
                self.runlog["iter "+str(self.iter+1)+" center diff"]=diff
                return True
        else:
            lose_center=True
            diff=np.nan
        self.runlog["iter "+str(self.iter+1)+" lose center"]=lose_center
        self.runlog["iter "+str(self.iter+1)+" center diff"]=diff    
        return False
    
    def GetFinalGroup_pre(self,list_of_orth):#parallelize!!!!!
        for i in list_of_orth:
            a=self.GetMinimalGroup(i)[1]
            centers=list(np.array(a)[0])
            cols=list(a.columns)
            for j in range(self.GetOrthGroupLength(i)):
                for k in range(len(cols)):
                    if j+1 in cols[k]:
                        self.finalgroup.append((i,j+1,centers[k]))

    def iterplus(self):
        self.iter+=1
    
    def GetFinalGroup(self):
        a=self.finalgroup[:]
        self.runlog["final group"]=pd.DataFrame({"orth_group":[i[0] for i in a],"num":[i[1] for i in a],"group":[i[2] for i in a]})
    
    def GetStartTime(self):
        self.starttime=datetime.datetime.now()
        self.runlog["program start time"]=self.starttime
    
    def GetTotalTime(self):
        self.runlog["total run time"]=str(datetime.datetime.now()-self.starttime)

    def Save_runlog(self,filepath):
        normal_dict = dict(self.runlog)
        filename=filepath+"/w="+str(self.w)+".pkl"
        with open(filename, 'wb') as f_obj:
            pickle.dump(normal_dict, f_obj)