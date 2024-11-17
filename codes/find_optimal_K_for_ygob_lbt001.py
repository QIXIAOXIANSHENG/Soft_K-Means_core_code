import numpy as np
import os
import copy
import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import matplotlib.pyplot as plt

df=pd.read_csv("D:/Grad study/lab/data/soft kmeans/gene selection/ygob_lbt001_all_4_species_with_orth_nogene.csv")
values = df.iloc[:,2:].T.values 
values = values.astype('float64') 
tool = MinMaxScaler(feature_range=(0, 1))
data = tool.fit_transform(values)
data=pd.DataFrame(data).T
data.insert(0, "orth_group", df["orth_group"])
data.insert(1, "num", df["num"])
data.columns=df.columns
sample_data_gene=df.columns[2:]
X=np.array(data[sample_data_gene])
distortions = []
inertias = []
mapping1 = {}
mapping2 = {}
K=range(1,20)
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