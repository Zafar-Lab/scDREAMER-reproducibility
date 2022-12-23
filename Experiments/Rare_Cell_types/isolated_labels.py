import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import silhouette_score
from sklearn.metrics import f1_score

# from scIB.clustering import opt_louvain
def opt_louvain(adata, label_key, cluster_key, function=None, resolutions=None,
                use_rep=None,
                inplace=True, plot=False, force=True, verbose=True, **kwargs):
    """
    params:
        label_key: name of column in adata.obs containing biological labels to be
            optimised against
        cluster_key: name of column to be added to adata.obs during clustering. 
            Will be overwritten if exists and `force=True`
        function: function that computes the cost to be optimised over. Must take as
            arguments (adata, group1, group2, **kwargs) and returns a number for maximising
        resolutions: list if resolutions to be optimised over. If `resolutions=None`,
            default resolutions of 20 values ranging between 0.1 and 2 will be used
        use_rep: key of embedding to use only if adata.uns['neighbors'] is not defined,
            otherwise will be ignored
    returns:
        res_max: resolution of maximum score
        score_max: maximum score
        score_all: `pd.DataFrame` containing all scores at resolutions. Can be used to plot the score profile.
        clustering: only if `inplace=False`, return cluster assignment as `pd.Series`
        plot: if `plot=True` plot the score profile over resolution
    """

    if verbose:
        print('Clustering...')

    if function is None:
        function = metrics.nmi

    if cluster_key in adata.obs.columns:
        if force:
            print(f"Warning: cluster key {cluster_key} already exists "
                  "in adata.obs and will be overwritten")
        else:
            raise ValueError(f"cluster key {cluster_key} already exists in " +
                             "adata, please remove the key or choose a different name." +
                             "If you want to force overwriting the key, specify `force=True`")

    if resolutions is None:
        n = 10
        resolutions = list(np.linspace(1, 2, num=n))

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []
    print('use rep:',use_rep)
    # if use_rep=='X_pca':
    #     sc.pp.neighbors(adata)
    # else:
    try:
        adata.uns['neighbors']
    except KeyError:
        if verbose:
            print('computing neighbours for opt_cluster')
        sc.pp.neighbors(adata, use_rep=use_rep)

    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key)
        noofclus = len(np.unique(adata.obs[cluster_key]))
        score = function(adata, label_key, cluster_key, **kwargs)
        if verbose:
            print(f'resolution: {res}, no of clus:{noofclus}  ,{function.__name__}: {score}')
        score_all.append(score)
        if score_max < score:
            score_max = score
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]

    if verbose:
        print(f'optimised clustering against {label_key}')
        print(f'optimal cluster resolution: {res_max}')
        print(f'optimal score: {score_max}')

    score_all = pd.DataFrame(zip(resolutions, score_all), columns=('resolution', 'score'))
    if plot:
        # score vs. resolution profile
        sns.lineplot(data= score_all, x='resolution', y='score').set_title('Optimal cluster resolution profile')
        plt.show()

    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, score_max, score_all
    else:
        return res_max, score_max, score_all, clustering

def silhouette(
        adata,
        group_key,
        embed,
        metric='euclidean',
        scale=True
):
    """
    wrapper for sklearn silhouette function values range from [-1, 1] with 1 being an ideal fit, 0 indicating
    overlapping clusters and -1 indicating misclassified cells

    :param group_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    :param scale: default True, scale between 0 (worst) and 1 (best)
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    asw = silhouette_score(
        X=adata.obsm[embed],
        labels=adata.obs[group_key],
        metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw



def isolated_labels(
        adata,
        label_key,
        batch_key,
        embed,
        isolated_labels,
        cluster=True,
        n=None,
        all_=False,
        verbose=True
):
    """
    score how well labels of isolated labels are distiguished in the dataset by
        1. clustering-based approach F1 score
        2. average-width silhouette score on isolated-vs-rest label assignment
    params:
        cluster: if True, use clustering approach, otherwise use silhouette score approach
        embed: key in adata.obsm used for silhouette score if cluster=False, or
            as representation for clustering (if neighbors missing in adata)
        n: max number of batches per label for label to be considered as isolated.
            if n is integer, consider labels that are present for n batches as isolated
            if n=None, consider minimum number of batches that labels are present in
        all_: return scores for all isolated labels instead of aggregated mean
    return:
        by default, mean of scores for each isolated label
        retrieve dictionary of scores for each label if `all_` is specified
    """

    scores = {}
    
    for label in isolated_labels:
        score = score_isolated_label(
            adata,
            label_key,
            label,
            embed,
            cluster,
            verbose=verbose
        )
        scores[label] = score

    if all_:
        return scores
    return np.mean(list(scores.values()))


def score_isolated_label(
        adata,
        label_key,
        label,
        embed,
        cluster=True,
        iso_label_key='iso_label',
        verbose=False
):
    """
    compute label score for a single label
    params:
        adata: anndata object
        label_key: key in adata.obs of isolated label type (usually cell label)
        label: value of specific isolated label e.g. cell type/identity annotation
        embed: embedding to be passed to opt_louvain, if adata.uns['neighbors'] is missing
        cluster: if True, compute clustering-based F1 score, otherwise compute
            silhouette score on grouping of isolated label vs all other remaining labels
        iso_label_key: name of key to use for cluster assignment for F1 score or
            isolated-vs-rest assignment for silhouette score
    """
    adata_tmp = adata.copy()

    def max_f1(adata, label_key, cluster_key, label, argmax=False):
        """cluster optimizing over largest F1 score of isolated label"""
        obs = adata.obs
        max_cluster = None
        max_f1 = 0
        for cluster in obs[cluster_key].unique():
            y_pred = obs[cluster_key] == cluster
            y_true = obs[label_key] == label
            f1 = f1_score(y_pred, y_true)
            if f1 > max_f1:
                max_f1 = f1
                max_cluster = cluster
        if argmax:
            return max_cluster
        return max_f1

    if cluster:
        # F1-score on clustering
        opt_louvain(
            adata_tmp,
            label_key,
            cluster_key=iso_label_key,
            label=label,
            use_rep=embed,
            function=max_f1,
            verbose=False,
            inplace=True
        )
        score = max_f1(adata_tmp, label_key, iso_label_key, label, argmax=False)
    else:
        # AWS score between label
        adata_tmp.obs[iso_label_key] = adata_tmp.obs[label_key] == label
        score = silhouette(adata_tmp, iso_label_key, embed)

    del adata_tmp

    if verbose:
        print(f"{label}: {score}")

    return score