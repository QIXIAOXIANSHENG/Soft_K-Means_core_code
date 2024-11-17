import pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
import numpy as np
import os
from scipy.spatial.distance import cdist
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import calinski_harabasz_score

def dunn_index(X, labels):
    """
    Calculate the Dunn Index for a clustering result.
    
    Parameters:
    X: ndarray
        Data points (samples x features)
    labels: array-like
        Cluster labels for each point in X
    
    Returns:
    float
        Dunn Index value
    """
    
    # Get the unique cluster labels
    unique_labels = np.unique(labels)
    
    # Initialize intra-cluster and inter-cluster distances
    intra_dists = []
    inter_dists = []

    # Calculate intra-cluster distances (within each cluster)
    for label in unique_labels:
        cluster_points = X[labels == label]
        if len(cluster_points) > 1:
            # Compute the pairwise distances within the cluster
            intra_dist = np.max(cdist(cluster_points, cluster_points))
            intra_dists.append(intra_dist)
    
    # Calculate inter-cluster distances (between different clusters)
    for i, label_a in enumerate(unique_labels[:-1]):
        cluster_a = X[labels == label_a]
        for label_b in unique_labels[i + 1:]:
            cluster_b = X[labels == label_b]
            # Compute the pairwise distances between cluster_a and cluster_b
            inter_dist = np.min(cdist(cluster_a, cluster_b))
            inter_dists.append(inter_dist)

    # Dunn Index: minimum inter-cluster distance / maximum intra-cluster distance
    return np.min(inter_dists) / np.max(intra_dists)



readdir=r"D:\Grad study\lab\data\soft kmeans\gene selection\ygob_libnorm_lbt005_cluster_result\plot_heatmap_csv file"
savedir=r'D:\Grad study\lab\data\soft kmeans\gene selection\ygob_libnorm_lbt005_cluster_result\scores'
reward=[0,0.05,0.1,0.2,0.5]
k_list=[3,4,5,7,9]
silhouette=[]
c_h=[]
d=[]
d_b=[]
for k in k_list:
    for w in reward:
        df=pd.read_csv(readdir+r"\k="+str(k)+r"\w="+str(w)+".csv")
        sil_score = silhouette_score(np.array(df.iloc[:,3:]), df["group"])
        ch_index = calinski_harabasz_score(np.array(df.iloc[:,3:]), df["group"])
        dunn=dunn_index(np.array(df.iloc[:,3:]),df["group"])
        db_index = davies_bouldin_score(np.array(df.iloc[:,3:]),df["group"])
        silhouette.append(sil_score)
        c_h.append(ch_index)
        d.append(dunn)
        d_b.append(db_index)


reward_str=[str(i) for i in reward]
res=pd.DataFrame({"k":list(np.repeat(k_list, len(reward))),"w":reward_str*len(k_list),"silhouette":silhouette,"C-H":c_h,"Dunn":d,"D-B":d_b})
res.to_csv(savedir+r'\stat_scores.csv', index=False) 