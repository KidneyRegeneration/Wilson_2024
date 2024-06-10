import os
import csv
import time

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd

from scipy.spatial.distance import cdist
from mpl_toolkits.mplot3d import Axes3D
from anndata import AnnData
from collections import defaultdict
from typing import List, Tuple
import matplotlib.backends.backend_pdf as pdf

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go

from pyvis.network import Network
import networkx as nx

def ordered_cluster_names(data: ad.AnnData) -> np.ndarray:
    """
    Returns the ordered cluster names from the 'cat_cluster' column of the AnnData object.
    
    Parameters:
        data (ad.AnnData): The AnnData object.
        
    Returns:
        np.ndarray: The ordered cluster names.
    """
    unique_values, indices = np.unique(data.obs['cat_cluster'], return_index=True)
    return unique_values[np.argsort(indices)]


def precompute_cluster_info(adata: ad.AnnData, clusters: List[str]) -> Tuple[List[np.ndarray], List[int]]:
    """
    Precomputes cluster indices and sizes for a given AnnData object and list of cluster names.
    
    Parameters:
        adata (ad.AnnData): The AnnData object containing the data.
        clusters (List[str]): The list of cluster names.
        
    Returns:
        Tuple[List[np.ndarray], List[int]]: A tuple containing the cluster indices and sizes.
    """
    cluster_indices = [np.where(adata.obs['cat_cluster'] == cluster)[0] for cluster in clusters]
    cluster_sizes = [len(indices) for indices in cluster_indices]
    return cluster_indices, cluster_sizes


def bootstrap_euclidean_interset(adata: ad.AnnData, reference_adata: ad.AnnData, n_iter: int) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Performs bootstrap sampling and computes Euclidean distances between cluster means of adata and reference_adata.
    
    Parameters:
        adata (ad.AnnData): The AnnData object containing the data.
        reference_adata (ad.AnnData): The reference AnnData object.
        n_iter (int): The number of iterations for bootstrap sampling.
        
    Returns:
        Tuple[np.ndarray, List[str], List[str]]: A tuple containing the distances, cluster names of adata, and cluster names of reference_adata.
    """
    clusters = ordered_cluster_names(adata)
    reference_clusters = ordered_cluster_names(reference_adata)
    
    cluster_indices, cluster_sizes = precompute_cluster_info(adata, clusters)
    cluster_indices_ref, cluster_sizes_ref = precompute_cluster_info(reference_adata, reference_clusters)
    
    n_vars = adata.n_vars
    n_clusts = len(clusters)
    n_clusts_ref = len(reference_clusters)
    
    distances = np.zeros((n_iter, n_clusts_ref, n_clusts))
    
    for i in range(n_iter):
        iter_means = np.zeros((n_clusts, n_vars))
        iter_means_ref = np.zeros((n_clusts_ref, n_vars))
        
        for c_number, (indices, size) in enumerate(zip(cluster_indices, cluster_sizes)):
            random_indices = np.random.choice(indices, size=size, replace=True)
            iter_means[c_number, :] = adata.X[random_indices, :].mean(axis=0)
        
        for c_number, (indices, size) in enumerate(zip(cluster_indices_ref, cluster_sizes_ref)):
            random_indices_ref = np.random.choice(indices, size=size, replace=True)
            iter_means_ref[c_number, :] = reference_adata.X[random_indices_ref, :].mean(axis=0)
        
        distances[i, :, :] = cdist(iter_means_ref, iter_means, metric='euclidean')
    
    return distances, clusters, reference_clusters


def bootstrap_euclidean_interset_orig(adata, reference_adata, n_iter):
    """
    Function calculates the Euclidean distance between the means of clusters
    in 'adata' and 'reference_adata' by bootstrap sampling, over 'n_iter' iterations.

    Args:
        adata: Anndata object containing the first dataset with cluster assignments.
        reference_adata: Anndata object containing the reference dataset with its 
                         own cluster assignments.
        n_iter (int): Number of bootstrap iterations to perform.

    Returns:
        numpy.ndarray: A 3D array of shape (n_iter, n_clusts_ref, n_clusts) containing 
                       the Euclidean distances calculated in each iteration.
        list: Names of the clusters in 'adata'.
        list: Names of the clusters in 'reference_adata'.
    """
    n_vars = adata.n_vars
    clusters = ordered_cluster_names(adata)
    n_clusts = len(clusters)
    reference_clusters = ordered_cluster_names(reference_adata)
    n_clusts_ref = len(reference_clusters)

    distances = np.zeros(( n_iter, n_clusts_ref, n_clusts ))

    for i in range(n_iter):
        
        iter_means = np.zeros((n_clusts, n_vars))
        iter_means_ref = np.zeros((n_clusts_ref, n_vars))

        cluster_indices = [np.where(adata.obs['cat_cluster'] == c)[0] for c in clusters]
        cluster_sizes = [len(indices) for indices in cluster_indices]
        random_indices = [np.random.choice(indices, size=size, replace=True) for indices, size in zip(cluster_indices, cluster_sizes)]
        
        for c_number, indices in enumerate(random_indices):
            iter_means[c_number, :] = adata.X[indices, :].mean(axis=0)
        
        cluster_indices_ref = [np.where(reference_adata.obs['cat_cluster'] == c)[0] for c in reference_clusters]
        cluster_sizes_ref = [len(indices) for indices in cluster_indices_ref]
        random_indices_ref = [np.random.choice(indices, size=size, replace=True) for indices, size in zip(cluster_indices_ref, cluster_sizes_ref)]
        
        for c_number, indices in enumerate(random_indices_ref): 
            iter_means_ref[c_number, :] = reference_adata.X[indices, :].mean(axis=0)

        distances[i, :, :] = cdist(iter_means_ref, iter_means, metric='euclidean')
    
    return distances, clusters, reference_clusters  


def filter_common_genes(adata1: ad.AnnData, adata2: ad.AnnData) -> Tuple[ad.AnnData, ad.AnnData]:
    common_genes = list(set(adata1.var_names) & set(adata2.var_names))
    if len(common_genes) == 0:
        raise ValueError("No common genes found between the datasets.")
    elif len(common_genes) < 50:
        raise ValueError("Number of common genes is less than 50.")
    return adata1[:, common_genes].copy(), adata2[:, common_genes].copy()


def make_hist(counts: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:     
    n_bins = int(len(counts) / 30)
    x_min = np.min(counts)
    x_max = np.max(counts)
    ys, bin_edges = np.histogram(counts, bins=n_bins, range=(x_min, x_max))
    if x_min == x_max: 
        raise ValueError("x_min is equal to x_max. Cannot create histogram.")
        
    bin_width = (x_max - x_min) / n_bins
    xs = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    return xs, ys / bin_width, (x_min + x_max) / 2.0


def remove_small_clusters(adata: ad.AnnData) -> Tuple[List[str], List[str]]:
    min_size: int = 50
    clusters: List[str] = ordered_cluster_names(adata)
    large_enough_clusters: List[str] = []
    removed_clusters: List[str] = []
    for cluster in clusters:
        num_obs: int = (adata.obs['cat_cluster'] == cluster).sum()
        print(f"{cluster} has {num_obs} observations")
        if num_obs > min_size:
            large_enough_clusters.append(cluster)
        else: 
            print(f"{cluster} Removed")
            removed_clusters.append(cluster)
    
    return large_enough_clusters, removed_clusters


def make_sankey_ordered(
    labels: List[str], 
    sources: List[str], 
    targets: List[str], 
    values: List[float], 
    title: str,
    output_file: str
) -> None:
    link_count = defaultdict(lambda: {'source': 0, 'target': 0})
    for src, tgt in zip(sources, targets):
        link_count[src]['source'] += 1
        link_count[tgt]['target'] += 1
        
    def sorting_criteria(link):
        src, tgt = link[0], link[1]
        src_count, tgt_count = link_count[src], link_count[tgt]
        is_SL = src_count['source'] == 1 and tgt_count['target'] == 1
        is_S_two_T_not_target = src_count['source'] == 2 and tgt_count['target'] == 0
        is_S_two_T_can_be_target = src_count['source'] == 2
        is_S_more_than_two_T_not_target = src_count['source'] > 2 and tgt_count['target'] == 0
        return (is_SL, is_S_two_T_not_target, is_S_two_T_can_be_target, is_S_more_than_two_T_not_target)
        
    sorted_links = sorted(zip(sources, targets, values), key=sorting_criteria, reverse=True)
     
    sorted_sources, sorted_targets, sorted_values = zip(*sorted_links)
    fig = go.Figure(go.Sankey(
        node={"label": labels, 'pad': 10, 'thickness': 15},
        link={"source": sorted_sources, "target": sorted_targets, "value": sorted_values}
    ))
    fig.update_layout(title_text=title, font_size=13, height=900)
    fig.write_html(output_file)


def generate_colors(n):
    """
    Generate a list of n distinct colors.
    
    Args:
    - n (int): The number of colors to generate.
    
    Returns:
    - list: A list of RGB color values.
    """
    colors = [mcolors.hsv_to_rgb((i/n, 1, 1)) for i in range(n)]
    return colors


def histograms_distances(
    distrs_d: np.ndarray, 
    clusters_x: np.ndarray, 
    clusters_y: np.ndarray, 
    output_file: str
) -> None:
    """
    Generate histograms for distances between clusters.

    Args:
        distrs_d (np.ndarray): Array of distances between clusters.
        clusters_x (np.ndarray): Array of labels for the x-axis clusters.
        clusters_y (np.ndarray): Array of labels for the y-axis clusters.
        output_file (str): Path to the output PDF file.

    Returns:
        None
    """
    rows_per_figure: int = 4
    cols_per_figure: int = 1
    colors: List[str] = generate_colors(len(clusters_x))
    pdf_pages: pdf.PdfPages = pdf.PdfPages(output_file)  # Create a PDF file to save the plots

    for fig_num in range(-(-len(clusters_y) // rows_per_figure)):
        fig, axs = plt.subplots(rows_per_figure, cols_per_figure, figsize=(20, 11), constrained_layout=True)

        start_num: int = fig_num * rows_per_figure
        end_num: int = min(start_num + rows_per_figure, len(clusters_y))

        for num_y in range(start_num, end_num):
            row: int = num_y - start_num
            axs[row].set_ylabel(clusters_y[num_y], fontsize=10)
            for num_x, clust_x in enumerate(clusters_x):
                if len(clusters_x) == np.shape(distrs_d)[1]:
                    xs, ys, x_annotate = make_hist(distrs_d[:, num_x, num_y])
                else:
                    raise ValueError("Length of clusters_y does not match the shape of distances.")

                axs[row].plot(xs, ys, linewidth=2, color=colors[num_x])
                axs[row].annotate(clusters_x[num_x], (x_annotate, max(ys)), color=colors[num_x], fontsize=13)

        # Hide unused subplots in the last figure
        for i in range(end_num - start_num, rows_per_figure):
            axs[i].axis('off')

        # Save the current figure to a file with fig_num appended to the filename
        output_file_with_fig_num = f"{output_file}_{fig_num}.pdf"
        fig.savefig(output_file_with_fig_num)

    plt.close('all')  # Close all figures

def qs_directed_interset(
    clusts_from: np.ndarray, 
    clusts_into: np.ndarray, 
    means: np.ndarray, 
    stds: np.ndarray
) -> pd.DataFrame:
    """
    Calculate the significance of alignments between clusters (q_ij).

    Args:
        clusts_from (np.ndarray): Array of clusters from which the distances are calculated.
        clusts_into (np.ndarray): Array of clusters into which the distances are calculated.
        means (np.ndarray): Array of mean distances between clusters.
        stds (np.ndarray): Array of standard deviations of distances between clusters.

    Returns:
        pd.DataFrame: DataFrame containing the calculated quality scores.
    """
    rows_list: List[Dict[str, Union[float, str]]] = []  
    clusts: np.ndarray = clusts_from.copy()
    for clust_number, clust in enumerate(clusts):
        clust_means: np.ndarray = means[:, clust_number].copy()
        clust_means = clust_means[np.nonzero(clust_means)] 
        clust_means = clust_means[np.argsort(clust_means)]
        coord_nn: int = np.where(means[:, clust_number] == clust_means[0])[0][0]
        
        rows_list.append({
            'q_ij': 0,
            '<d_ij>': 0,
            'std_<d_ij>': 0,
            'oc': clust,
            'c_j': clusts_into[coord_nn],
        })

        nnn: int = 0
        clust_means = clust_means[1:]
        
        while nnn < len(clust_means):    
            coord_nnn: int = np.where(means[:, clust_number] == clust_means[nnn])[0][0]
            dij: float = means[coord_nnn, clust_number] - means[coord_nn, clust_number]
            stdij: float = np.sqrt(stds[coord_nn, clust_number]**2 + (stds[coord_nnn, clust_number]**2))
            
            rows_list.append({
                'q_ij': dij / stdij,
                '<d_ij>': dij,
                'std_<d_ij>': stdij,
                'oc': clust,
                'c_j': clusts_into[coord_nnn],
            })

            nnn += 1

    return pd.DataFrame(rows_list)

def absolute_distances_clusts(
    clusts: np.array, 
    clusters_ref: np.ndarray, 
    distances: np.array, 
    file_path: str
) -> None:
    """
    Calculate the absolute distances between clusters.

    Args:
        clusts (np.array): Array of clusters from which the distances are calculated.
        clusters_ref (np.ndarray): Array of reference clusters, into which the distances are calculated.
        distances (np.array): Array of distances between clusters.
        file_path (str): Path to the output file.

    Returns:
        None
    """
    means = np.mean(distances, axis=0)    
    stds = np.std(distances, axis=0) 
    rows_list = []  # List to store rows of the DataFrame
    for clust_number, clust in enumerate(clusts):
        clust_means = means[:, clust_number]
        clust_means = clust_means[np.nonzero(clust_means)] 
        clust_means = clust_means[np.argsort(clust_means)]
        coord_nn = np.where(means[:, clust_number] == clust_means[0])[0][0]
        
        rows_list.append({
            'q_ij': 0,
            '<d_ij>': means[coord_nn, clust_number],
            'oc': clust,
            'c_j': clusters_ref[coord_nn]
        })

        nnn = 0
        clust_means = clust_means[1:]
        
        while nnn < len(clust_means):    
            coord_nnn = np.where(means[:, clust_number] == clust_means[nnn])[0][0]
            dij = means[coord_nnn, clust_number] - means[coord_nn, clust_number]
            stdij = np.sqrt(stds[coord_nn, clust_number]**2 + (stds[coord_nnn, clust_number]**2))
            
            rows_list.append({
                'q_ij': dij / stdij,
                '<d_ij>': means[coord_nnn, clust_number],
                'oc': clust,
                'c_j': clusters_ref[coord_nnn]
            })

            nnn += 1

    df = pd.DataFrame(rows_list)
    df.to_csv(file_path, index=False)

def network_to_html(labels: np.ndarray[str], sources: np.ndarray, targets: np.ndarray, output_file: str) -> None:
    G = nx.DiGraph()
    for source, target in zip(sources, targets):
        G.add_edge(labels[source], labels[target])
    net = Network(height="750px", width="100%", directed=False, notebook=False)
    net.from_nx(G)
    net.write_html(output_file)

def alignment_table_interset_oneway(cut_off: float, dists: np.ndarray, clusters: np.ndarray, clusters_ref: np.ndarray, output_file: str) -> None: 
    qs_df = qs_directed_interset(clusters, clusters_ref, np.mean(dists, axis=0), np.std(dists, axis=0))
    df_sorted = qs_df.sort_values(by=['q_ij'], ascending=True, axis=0)

    aligned_rows = df_sorted[df_sorted['q_ij'] <= cut_off]
    aligned_indices = aligned_rows.index
    combined_rows = df_sorted.loc[aligned_indices]
    
    with open(output_file, 'w') as f:
        for index, clust in enumerate(clusters):
            clust_select_rows = combined_rows[combined_rows['oc'] == clust]
            
            # Skip clusters with no selected rows
            if clust_select_rows.empty:
                continue
            
            # Write cluster identifier
            f.write(f"Cluster: {clust}\n")
            f.write(clust_select_rows.drop(columns=['oc']).to_string(index=False))
            f.write("\n\n")  # Add some space before the next cluster for readability   
    
def alignment_diagram_interset_oneway(
    cut_off: float, 
    distances: np.ndarray, 
    clusters: np.ndarray, 
    reference_clusters: np.ndarray, 
    remove_this: np.ndarray, 
    output_file_sankey: str,
    output_file_network: str
) -> None:
    """
    Generate sankey and network diagrams of alignmnets between clusters in 2 datasets.

    Args:
        cut_off (float): significance cut-off for alignment
        distances (np.ndarray): Array of distances between clusters
        clusters (np.ndarray): Array of clusters from which the distances were calculated
        reference_clusters (np.ndarray): Array of reference clusters, into which the distances were calculated
        remove_this (np.ndarray): Array of small clusters to remove.
        output_file_sankey (str): Output file path for the Sankey diagram.
        output_file_network (str): Output file path for the network diagram.

    Returns:
        None
    """
    means = np.mean(distances, axis=0)  
    
    targets = []
    sources = []
    values = []
    
    qs_df = qs_directed_interset(clusters.copy(), 
                                reference_clusters.copy(), 
                                np.mean(distances, axis=0), 
                                np.std(distances, axis=0))

    for c, cluster in enumerate(clusters):        
        
        matching_rows = qs_df[qs_df['oc'] == cluster]        
        filtered_rows = matching_rows[matching_rows['q_ij'] <= cut_off]
        
        for _, row in filtered_rows.iterrows():
            aligned_cluster_index = np.where(reference_clusters==row['c_j'])[0][0]
            targets.append(len(clusters) + aligned_cluster_index)
            sources.append(c)
            values.append(1/means[aligned_cluster_index, c])

    labels = [*clusters , *reference_clusters]

    unique_targets = np.unique(targets)
    unique_sources = np.unique(sources)
    
    x_list = np.append(0.1 * np.ones(len(unique_sources)), 0.6 * np.ones(len(unique_targets))).tolist()
    y_list = np.round(np.append(np.linspace(0.01, 0.9, len(unique_sources)), np.linspace(0.01, 0.9, len(unique_targets))), decimals=3).tolist()

    if len(remove_this) == 0:  
        title = "CAT. "
    else: 
        title = f"list of clusters removed (size<=50) is in the corresponding txt file. "            
                
    make_sankey_ordered(labels, sources, targets, values, title, output_file_sankey)
    network_to_html(labels, sources, targets, output_file_network)

def create_file_paths(type_test, test_str, base_dir_name="CAT_test_results"):
    # Moving one level up from the current working directory
    parent_dir = os.path.join(os.getcwd(), os.pardir)
    parent_dir = os.path.abspath(parent_dir)  # Normalize the path
    
    base_dir = os.path.join(parent_dir, base_dir_name)
    test_dir = os.path.join(base_dir, f"{test_str}")
    # Ensure the directories exist
    os.makedirs(test_dir, exist_ok=True)
    
    file_paths = {
        'removed_clusters_file': os.path.join(test_dir, "_removed_clusters.txt"),
        'histograms_file': os.path.join(test_dir, f"{type_test}_histograms"),
        'alignment_table_file': os.path.join(test_dir, f"{type_test}_alignment_table.csv"),
        'alignment_diagram_file': os.path.join(test_dir, f"{type_test}_sankey.html"),
        'absolute_distances_file': os.path.join(test_dir, f"{type_test}_absolute_distances.csv"),
        'network_diagram_file': os.path.join(test_dir, f"{type_test}_network.html")
    }
    return file_paths


def save_removed_clusters_text(file_path: str, removed_clusters: list):
    unique_removed_clusters = set(removed_clusters)  # Using set to get unique values
    with open(file_path, 'w') as f:
        for cluster in sorted(unique_removed_clusters):  # Optional sorting for consistency
            f.write(f"{cluster}\n")


def process_clusters(adata):
    clusters, removed_clusters = remove_small_clusters(adata)
    mask = adata.obs['cat_cluster'].isin(clusters)
    adata = adata[mask].copy()
    return adata, removed_clusters


def CAT_one_vs_all(adata_dict, n_iter, test_str, cut_off): 
    for key, adata_one in adata_dict.items():
        print(key)
        adata_list_excluding_current = [adata for k, adata in adata_dict.items() if k != key]
        adata_all = sc.concat(adata_list_excluding_current, join='inner')
        test_type = f"{key}_vs_three_others"
        adata_one, adata_all = filter_common_genes(adata_one, adata_all)
        adata_one, removed_clusters_one = process_clusters(adata_one)
        adata_all, removed_clusters_all = process_clusters(adata_all)
        assert np.all(adata_one.var_names == adata_all.var_names), "Gene intersection between two datasets don't match!"
        distances, cl, cl_ref = bootstrap_euclidean_interset(adata_one, adata_all, n_iter)

        file_paths = create_file_paths(test_type, test_str, base_dir_name="CAT_test_results")
        combined_removed_clusters = removed_clusters_one + removed_clusters_all
        save_removed_clusters_text(file_paths['removed_clusters_file'], combined_removed_clusters)
        
        histograms_distances(distances.copy(), cl_ref.copy(), cl.copy(), file_paths['histograms_file'])

        alignment_table_interset_oneway(cut_off, distances.copy(), cl.copy(), cl_ref.copy(), file_paths['alignment_table_file'])
        
        absolute_distances_clusts(cl.copy(), cl_ref.copy(), distances.copy(), file_paths['absolute_distances_file'])

        alignment_diagram_interset_oneway(cut_off, 
        distances.copy(), 
        cl.copy(), 
        cl_ref.copy(), 
        combined_removed_clusters, 
        file_paths['alignment_diagram_file'], 
        file_paths['network_diagram_file'])


def CAT_consecutive(adata_dict, n_iter, test_str, cut_off):
    # Convert dictionary keys and values to lists to ensure we can access them by index
    keys_list = list(adata_dict.keys())
    adata_list = list(adata_dict.values())
    
    for i in range(len(adata_list) - 1):
        # Current and next adata object and their keys
        key_one, adata_one = keys_list[i], adata_list[i]
        key_two, adata_two = keys_list[i + 1], adata_list[i + 1]
        
        # Filter for common genes between the two consecutive adata objects
        adata_one, adata_two = filter_common_genes(adata_one, adata_two)
        
        # Process clusters for both adata objects
        adata_one, removed_clusters_one = process_clusters(adata_one)
        adata_two, removed_clusters_two = process_clusters(adata_two)
        combined_removed_clusters = removed_clusters_one + removed_clusters_two
        # Ensure that both adata objects have the same variable (gene) names after processing
        assert np.all(adata_one.var_names == adata_two.var_names), "Gene intersection between two datasets don't match!"
        
        # Calculate distances between the two datasets
        distances, cl, cl_ref = bootstrap_euclidean_interset(adata_one, adata_two, n_iter)
        def save_results(test_type, distances, cl, cl_ref, combined_removed_clusters):
            file_paths = create_file_paths(test_type, test_str, base_dir_name="CAT_test_results")
            histograms_distances(distances.copy(), cl_ref.copy(), cl.copy(), file_paths['histograms_file'])
            alignment_table_interset_oneway(cut_off, distances.copy(), cl.copy(), cl_ref.copy(), file_paths['alignment_table_file'])
            absolute_distances_clusts(cl.copy(), cl_ref.copy(), distances.copy(), file_paths['absolute_distances_file'])
            alignment_diagram_interset_oneway(cut_off, 
            distances.copy(), 
            cl.copy(), 
            cl_ref.copy(), 
            combined_removed_clusters, 
            file_paths['alignment_diagram_file'], 
            file_paths['network_diagram_file'])
        
        save_results(f"from_{key_one}_to_{key_two}", distances, cl, cl_ref, combined_removed_clusters)
        transposed_distances = np.transpose(distances, (0, 2, 1))

        save_results(f"from_{key_two}_to_{key_one}", transposed_distances, cl_ref, cl, combined_removed_clusters)


def main_bootstrap(filenames: List[str], obs_names: List[str], n_iter: int, cut_off: float, test_structure: dict[str, dict[str, List[str]]] ) -> None:
    
    tests = list(test_structure.keys())
    stages = ["Stage_1", "Stage_2", "Stage_3", "Stage_4"]  # Corrected

    for test in tests:
        print(f"Processing test: {test}")
        start_time = time.time()

        adata_dict = {}  # Dictionary to store each adata object

        for file_name, obs_name, stage_name, condition_name in zip(filenames, obs_names, stages, (test_structure[test]["conditions"])):
            adata = ad.read_h5ad(file_name)     
            adata = adata[adata.obs['condition'] == condition_name].copy()
            adata.obs['cat_cluster'] = adata.obs[obs_name].astype(str) + '_' + adata.obs['Stage'].astype(str) + '_' + adata.obs['condition'].astype(str)
            adata_dict[stage_name] = adata
        
        CAT_one_vs_all(adata_dict, n_iter, test, cut_off)
        CAT_consecutive(adata_dict, n_iter, test, cut_off)
    elapsed_time = time.time() - start_time
    print(f"Processing for {test} completed in {elapsed_time:.2f} seconds")
            

test_structure = {
    "Test1": {"conditions": ["3c", "3c4F", "3c4F", "3c4F"]},
    "Test2": {"conditions": ["3c", "3c4F", "3c4FG11", "3c4FG11"]},
    "Test3": {"conditions": ["3c", "3c4F200", "3c4F200", "3c4F200"]},
    "Test4": {"conditions": ["3c", "3c4FA", "3c4FA", "3c4FA"]},
    "Test5": {"conditions": ["3c", "3c4FC", "3c4FC", "3c4FC"]},
    "Test6": {"conditions": ["3c", "3c4FL", "3c4FL", "3c4FL"]},
    "Test7": {"conditions": ["3c", "3c4E", "3c4E", "3c4E"]},
    "Test8": {"conditions": ["4c", "4c4F", "4c4F", "4c4F"]},
    "Test9": {"conditions": ["4c", "4c4F", "4c4FG11", "4c4FG11"]},
    "Test10": {"conditions": ["4c", "4c4F200", "4c4F200", "4c4F200"]},
    "Test11": {"conditions": ["4c", "4c4FA", "4c4FA", "4c4FA"]},
    "Test12": {"conditions": ["4c", "4c4FC", "4c4FC", "4c4FC"]},
    "Test14": {"conditions": ["4c", "4c4E", "4c4E", "4c4E"]},
    "Test15": {"conditions": ["5c", "5c4F", "5c4F", "5c4F"]},
    "Test16": {"conditions": ["5c", "5c4F", "5c4FG11", "5c4FG11"]},
    "Test17": {"conditions": ["5c", "5c4F200", "5c4F200", "5c4F200"]},
    "Test18": {"conditions": ["5c", "5c4FA", "5c4FA", "5c4FA"]},
    "Test19": {"conditions": ["5c", "5c4FC", "5c4FC", "5c4FC"]},
    "Test20": {"conditions": ["5c", "5c4FL", "5c4FL", "5c4FL"]},
    "Test21": {"conditions": ["5c", "5c4E", "5c4E", "5c4E"]},
}


n_iter = 1500
cut_off = 1.6

#distance_metric = "Eucledian"

filenames = ["/Users/alisa/myenv/data/stage1_240124.h5ad", "/Users/alisa/myenv/data/stage2_240124.h5ad", "/Users/alisa/myenv/data/stage3_240124.h5ad", "/Users/alisa/myenv/data/stage4_240124.h5ad"]
obs_names = ['Ann', 'ann', 'ann', 'ann']

print("Current Working Directory:", os.getcwd())

main_bootstrap(filenames, obs_names, n_iter, cut_off, test_structure)