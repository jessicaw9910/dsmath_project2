#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import kneighbors_graph
import networkx as nx
from itertools import count
from scipy import sparse
from scipy import linalg

def create_distance(mx_input):
    '''Create adjacency matrix from input
    
    Params:
    - mx_input (numpy.ndarray): centered data matrix

    Output:
    - mx_distance (numpy.ndarray): distance matrix
    '''
    mx_distance = np.linalg.norm(mx_input[:, None, :] - mx_input[None, :, :], axis=-1)
    return mx_distance

def plot_epsilon(mx_input, eps_start=10, eps_end=40, eps_step=5, n_cols=3, 
                 title='plots/epsilon', save=True, img_size=(10,5)):
    '''
    Params:
    - mx_input (numpy.ndarray): centered data matrix
    - eps_start (int): np.arange start value of epsilon
    - eps_end (int): np.arange end value of epsilon (not included)
    - eps_step (int): np.arange step value of epsilon
    - n_cols (int): number of columns in figure - always two rows of plots
    - title (str): output name, including folder name
    - save (bool): saves figure if True
    - fig_size (tuple of floats): output figure size
    '''
    epsilon_list = list(np.arange(eps_start, eps_end, eps_step))
    
    fig, ax = plt.subplots(2, n_cols, figsize=img_size, tight_layout=True)
    
    for i, epsilon in zip(range(len(epsilon_list)), epsilon_list):
        ## first row of figure
        if epsilon <= epsilon_list[n_cols-1]:
            ax[0, i].hist(np.sum(np.where((mx_input <= epsilon) & (mx_input != 0), 1, 0), axis=1));
            ax[0, i].set_title('Epsilon: ' + str(epsilon));
        ## second row of figure
        else:
            ax[1, i-n_cols].hist(np.sum(np.where((mx_input <= epsilon) & (mx_input != 0), 1, 0), axis=1));
            ax[1, i-n_cols].set_title('Epsilon: ' + str(epsilon));
    
    fig.add_subplot(111, frame_on=False);
    plt.tick_params(labelcolor="none", bottom=False, left=False);
    plt.xlabel('Number of neighbors');
    plt.ylabel('Number of samples');
    
    if save:
        plt.savefig(title + '.svg',bbox_inches='tight');
    plt.show();    

def create_adjacency(mx_input, knn=True, n=10, epsilon=10.0):
    '''Create adjacency matrix from input

    Params:
    - mx_input (numpy.ndarray): centered data matrix
    - method (bool): 'knn' method if true else 'epsilon'
    - n (int): number of nearest neighbors included if knn=True
    - epsilon (float): distance threshold if knn=False

    Output:
    - mx_adjacency (numpy.ndarray): adjacency matrix
    '''
    if knn:
        mx_connectivity = kneighbors_graph(X=mx_input, n_neighbors=n, mode='connectivity')
        mx_adjacency = (1/2)*(mx_connectivity + mx_connectivity.T)
        mx_adjacency = mx_connectivity.todense()
    else:
        mx_distance = np.linalg.norm(mx_input[:, None, :] - mx_input[None, :, :], axis=-1)
        mx_adjacency = np.where((mx_distance <= epsilon) & (mx_distance != 0), 1, 0)
        
    return mx_adjacency

def create_graph(mx_adjacency, list_names):
    '''Create networkx graph object

    Params:
    - mx_adjacency (numpy.ndarray): adjacency matrix
    - list_names (list): list of items to rename nodes with

    Output:
    - G (networkx.classes.multigraph.MultiGraph): graph object corresponding to adjacency matrix
    '''
    G = nx.from_numpy_matrix(mx_adjacency, create_using=nx.MultiGraph)
    mapping = dict(zip(list(G.nodes()), list_names))
    G = nx.relabel_nodes(G, mapping)

    return G

def add_annotation(graph, df_annotate, col_name, col_annotate, annot_name):
    ''' Add annotations from drug/cell dataframes to graph object

    Params:
    - graph (networkx.classes.multigraph.MultiGraph): relevant graph
    - df_annotate (pd.Dataframe): dataframe containing relevant annotations
    - col_name (str): column name that corresponds to node name
    - col_annotate (str): column name that corresponds to info being annotated
    - annot_name (str): name of annotation added to graph object

    Output:
    - graph (networkx.classes.multigraph.MultiGraph): modified graph object
    '''
    nodes = graph.nodes()
    temp_list = []

    for node in list(nodes):
        temp_list.append(df_annotate.loc[df_annotate[col_name] == node, col_annotate].iloc[0])

    temp_dict = dict(zip(list(nodes), temp_list))
    nx.set_node_attributes(graph, temp_dict, annot_name)

    return graph


def return_graph(mx_adj, list_node, df_info, list_colname, list_colannot, list_annot):
    ''' Combine graph create and annotation formulas
    Params:
    - mx_adj (numpy.ndarray): adjacency matrix
    - list_node (list): list of items to rename nodes with
    - df_info (pd.DataFrame): dataframe containing relevant annotations
    - list_colname (list of strings): list of column names that corresponds to node name
    - list_colannot (list of strings): list of column names that corresponds to info being annotated
    - list_annot (list of strings): list of names of annotations added to graph object

    Output:
    - graph (networkx.classes.multigraph.MultiGraph): relevant graph with annotations

    '''
    graph = create_graph(mx_adj, list_node)
    
    for col_name, col_annot, annot in zip(list_colname, list_colannot, list_annot):
        graph = add_annotation(graph, df_info, col_name, col_annot, annot)
        
    return graph

## code from https://stackoverflow.com/questions/28910766/python-networkx-set-node-color-automatically-based-on-number-of-attribute-opt
def plot_graph(graph, attr, title='plots/graph', remove=False, save=True, img_size=(40, 20)):
    '''Return plot of graph

    Params:
    - graph (networkx.classes.multigraph.MultiGraph): relevant graph
    - attr (str): 
    - title (str): output name, including folder name
    - save (bool): saves figure if True
    - fig_size (tuple of floats): output figure size

    '''
    ## remove unconnected nodes
    if remove:
        graph.remove_nodes_from(list(nx.isolates(graph)))

    ## color graph by annotation
    nodes = graph.nodes()
    groups = set(nx.get_node_attributes(graph, attr).values())
    mapping = dict(zip(sorted(groups), count()))
    colors = [mapping[graph.node[n][attr]] for n in nodes]

    fig, ax = plt.subplots(figsize=img_size)
    pos = nx.spring_layout(graph)
    ec = nx.draw_networkx_edges(graph, pos, alpha=0.2);
    nc = nx.draw_networkx_nodes(graph, pos, nodelist=nodes, node_color=colors, cmap=plt.cm.jet);
    nx.draw_networkx_labels(graph, pos);
    plt.axis('off');

    ## create legend
    bounds = np.arange(0, len(groups)+1, 1)
    loc = bounds + 0.5
    cbar = plt.colorbar(nc, spacing='proportional', ticks=bounds, boundaries=bounds);
    cbar.set_ticks(loc);
    cbar.set_ticklabels(sorted(list(groups)));

    if save:
        plt.savefig(title + '.svg',bbox_inches='tight');
    plt.show();

def create_laplacian(mx_adjacency):
    ''' Create graph Laplacian
    
    Params:
    mx_adjacency (numpy.ndarray): adjacency matrix

    Output:
    mx_laplacian (numpy.ndarray): graph laplacian matrix
    '''

    mx_laplacian = sparse.csgraph.laplacian(csgraph=mx_adjacency, normed=True)
    return mx_laplacian

def compute_spectrum(mx_laplacian):
    '''Compute eigenvalues and eigenvectors and project onto the real numbers.
    
    Params:
    mx_laplacian (numpy.ndarray): graph laplacian matrix
    
    Output:
    eigenvals_sorted
    eigenvcts_sorted
    '''
    eigenvals, eigenvcts = linalg.eig(mx_laplacian)
    eigenvals = np.real(eigenvals)
    eigenvcts = np.real(eigenvcts)

    idx = np.argsort(eigenvals)
    eigenvals_sorted = eigenvals[idx]
    eigenvcts_sorted = eigenvcts[:, idx]

    return eigenvals_sorted, eigenvcts_sorted