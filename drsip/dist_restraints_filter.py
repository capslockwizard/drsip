"""
=========================================================
Distance Restraints Filter (:mod:`drsip.dist_restraints`)
=========================================================

Module contains the function implementing the distance restraints
filter.

Functions
---------

.. autofunction:: dist_restraints_filter
"""
from scipy import stats
from MDAnalysis.analysis import distances


def dist_restraints_filter(static_coord, mobile_coord, dist_mat_idx, dist_rest_data, cutoff=0.3):
    """Distance restraints filter.

    Pearson and Spearman correlation between the given distance
    restraints and the distances between C-alpha atoms of the
    corresponding residue pairs in docking poses are used to filter the
    docking poses.

    Parameters
    ----------
    static_ori_coord, mobile_ori_coord : np.array
        Nx3 coordinate array of the static monomer and mobile monomer,
        respectively, from the docking pose. Where N is the number of
        atoms.
    dist_mat_idx : list
        Contains 2 lists. The lists represents the row (first list) and
        column (second) indices (0-indexed) of the distance matrix that
        corresponds to the residue pairs in the distance restraints
        table.
    dist_rest_data : np.array
        Array containing the distance restraints.
    cutoff : float, optional
        The distance restraints cutoff. If the docking pose's Spearman
        and Pearson correlation <= cuttoff, the docking pose fails to
        pass the filter.

    Returns
    -------
    list
        Contains the distance between the C-alpha atoms (np.array),
        Spearman correlation (float), Pearson correlation (float) and
        pass status (bool). The pass status is True when the docking
        pose passes the filter and vice versa.
    """
    dist_mat = distances.distance_array(static_coord, mobile_coord)

    current_distances = dist_mat[tuple(dist_mat_idx)]
    spearman_correl = stats.pearsonr(stats.rankdata(
        dist_rest_data), stats.rankdata(current_distances))[0]
    pearson_correl = stats.pearsonr(dist_rest_data, current_distances)[0]
    pass_status = (spearman_correl > cutoff) and (pearson_correl > cutoff)

    return [current_distances, spearman_correl, pearson_correl, pass_status]


def dist_restraints_filter_pearson_only(static_coord, mobile_coord, dist_mat_idx, dist_rest_data, cutoff=0.3):
    dist_mat = distances.distance_array(static_coord, mobile_coord)

    current_distances = dist_mat[tuple(dist_mat_idx)]
    pearson_correl = stats.pearsonr(dist_rest_data, current_distances)[0]
    pass_status = pearson_correl > cutoff

    return [current_distances, pearson_correl, pass_status]
