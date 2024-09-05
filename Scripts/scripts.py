import numpy as np
import pandas as pd


###calculate markers' scores from a dictionary
def marker_score(markers_dict, adata, N_samples=100, random_seed=42):
    np.random.seed(random_seed)
    markers_list = []
    N_genes = adata.shape[1]
    random_genes = np.unique( np.random.randint(low=0, high=N_genes, size=N_samples) )
    gene_names = adata.var_names[random_genes]
    for i in markers_dict:
        markers_list.append(f'{i}_score')
        adata.obs[f'{i}_score'] = np.array(np.mean(adata[:,markers_dict[i]].X,1) - np.mean(adata[:,gene_names].X,(0,1)))
    return markers_list, adata

###function to rename clusters from a dictionary
def rename_clusters(names_dict, names_obs):
    clusters = pd.Categorical(names_obs)
    clusters=clusters.rename_categories(names_dict)
    cluster_array = np.array(clusters)
    split_array = [ i.split('.')[0] for i in cluster_array ]
    clusters = pd.Categorical(split_array)
    return clusters

###use scores instead of manual names (as above) to rename clusters
def clustersByScores(adata, markers_scores, leidenClusters):
    clusters = pd.Categorical(leidenClusters)
    scoresTable = adata.obs[markers_scores]
    clusterUnique = np.unique(leidenClusters)
    newNames = pd.Series(index=leidenClusters)
    for CLST in clusterUnique:
        meanScores = np.mean( scoresTable.loc[leidenClusters==CLST,:], 0)
        newId = meanScores.index[ np.argmax(meanScores) ].split('_')[0]
        newNames[CLST] = newId
    return(pd.Categorical(newNames))