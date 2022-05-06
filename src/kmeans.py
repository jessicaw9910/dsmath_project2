#!/usr/bin/env python3

import numpy as np
from sklearn.cluster import KMeans
from kneed import KneeLocator
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

def find_kmeans(mx_input, dict_kwarg, n_clust=30):
    sse = []
    silhouette = []

    for k in range(1, n_clust+1):
        kmeans = KMeans(n_clusters=k, **dict_kwarg)
        kmeans.fit(mx_input)
        sse.append(kmeans.inertia_)
        if k > 1:
            score = silhouette_score(mx_input, kmeans.labels_)
            silhouette.append(score)

    kl = KneeLocator(range(1, n_clust+1), sse, curve="convex", direction="decreasing")
    n_clust = kl.elbow

    kmeans = KMeans(n_clusters=n_clust, **dict_kwarg)
    kmeans.fit(mx_input)
    
    return kmeans, sse, silhouette

def plot_kmeans(kmeans_cell, kmeans_drug, sse_cell, sse_drug,
                title='knn', folder='plots/', save=True, img_size=(12,4)):
    fig, ax = plt.subplots(1, 2, figsize=img_size);
    ## cell
    plt.subplot(1, 2, 1);
    plt.plot(range(1, 31), sse_cell);
    plt.xlabel('Number of clusters');
    plt.ylabel('SSE');
    plt.vlines(x=len(np.unique(kmeans_cell.labels_)), ymin=0, ymax=max(sse_cell), colors='r', linestyle='dashed');
    plt.title('Cell');
    ## drug
    plt.subplot(1, 2, 2);
    plt.plot(range(1, 31), sse_drug);
    plt.xlabel('Number of clusters');
    plt.ylabel('SSE');
    plt.vlines(x=len(np.unique(kmeans_drug.labels_)), ymin=0, ymax=max(sse_drug), colors='r', linestyle='dashed');
    plt.title('Drug');

    if save:
        plt.savefig(folder + title + '.svg',bbox_inches='tight');
    plt.show();