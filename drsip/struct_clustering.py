"""
=====================================================
Structure Clustering (:mod:`drsip.struct_clustering`)
=====================================================

Module contains functions implementing the structure clustering part of
the protocol.

Functions
---------

.. autofunction:: min_dist
.. autofunction:: calc_MSD
.. autofunction:: calc_MSD_mob_static
.. autofunction:: cluster_poses
.. autofunction:: calc_rot_static_coord
.. autofunction:: cal_rmsd_between_poses_membrane
.. autofunction:: cal_rmsd_between_poses_soluble
"""

from numba import jit, float32
import scipy as sp
import scipy.cluster
import scipy.spatial.distance
import numpy as np
import drsip_common


def min_dist(dist_mat_1, dist_mat_2):
    """Compare the elements of 2 matrices and keep the minimum values.

    Values are compared element-wise and the minimum of the 2 values is
    returned.

    Parameters
    ----------
    dist_mat_1, dist_mat_2 : np.array
        NxN distance matrices between the 2 monomers in the docking
        pose. Where N are the number of atoms.

    Returns
    -------
    np.array
        Returns a new distance matrix containing the minimum values for
        each element.
    """
    temp_dist_mat = dist_mat_1.copy()
    swap_idx = np.where(dist_mat_1 > dist_mat_2)
    temp_dist_mat[swap_idx] = dist_mat_2[swap_idx]

    return temp_dist_mat


@jit(float32[:](float32[:, :], float32[:, :, :]), nopython=True)
def calc_MSD(ref_coord, mobile_coords):
    """Compute the mean squared deviation (MSD).

    MSDs are computed between the reference and the mobile monomers
    from different poses.

    Parameters
    ----------
    ref_coord : np.array
        Nx3 coordinate matrix of the reference monomer with N atoms.
    mobile_coords : np.array
        MxNx3 coordinate matrices of mobile monomers with N atoms
        from M poses.

    Returns
    -------
    np.array
        M-dim vector containing the MSDs between the static and mobile
        monomers.
    """
    num_coords = mobile_coords.shape[0]
    num_atoms = mobile_coords.shape[1]
    MSDs = np.zeros(num_coords, dtype=np.float32)

    for coord_idx in xrange(num_coords):
        MSDs[coord_idx] = np.sum(
            (mobile_coords[coord_idx] - ref_coord)**2)

    return MSDs/num_atoms


@jit(float32[:, :](float32[:, :, :], float32[:, :, :]), nopython=True)
def calc_MSD_mob_static(rotated_static_coords, mobile_coords):
    """
    Compute the MSD between mobile and static monomers.

    Parameters
    ----------
    rotated_static_coords : np.array
        MxNx3 coordinate matrix of the static monomers with N atoms
        from the M poses after superimposing the mobile monomer to the
        static monomer.
    mobile_coords : np.array
        MxNx3 coordinate matrices of the mobile monomers with N atoms
        from M poses.

    Returns
    -------
    np.array
        MxM matrix of the MSDs between all pairs of poses.
    """
    num_poses = rotated_static_coords.shape[0]
    MSDs = np.zeros((num_poses, num_poses), dtype=np.float32)

    # Compute the first and last rows of the MSD
    MSDs[0, 1:] = calc_MSD(rotated_static_coords[0], mobile_coords[1:])
    MSDs[num_poses-1, 0:num_poses -
         1] = calc_MSD(rotated_static_coords[num_poses-1],
                       mobile_coords[:num_poses-1])

    # For the other rows:
    for row_idx in xrange(1, num_poses-1):
        # To avoid comparing the same pose (diagonal), split the
        # calculation into two parts:

        # Lower half of the MSD
        MSDs[row_idx, 0:row_idx] = calc_MSD(
            rotated_static_coords[row_idx], mobile_coords[:row_idx])

        # Upper half of the MSDs
        MSDs[row_idx, row_idx +
             1:] = calc_MSD(rotated_static_coords[row_idx], mobile_coords[row_idx+1:])

    return MSDs


def cluster_poses(dist_mat, cutoff=12.0, criterion='distance'):
    """Cluster docking poses

    Docking poses are clustered with hierarchical clustering and
    average linkage criteria. The distance matrix contains the RMSD
    between the poses. See
    :py:meth:`cal_rmsd_between_poses_membrane <drsip.struct_clustering.cal_rmsd_between_poses_membrane>` for HoTPs
    or :py:meth:`cal_rmsd_between_poses_soluble <drsip.struct_clustering.cal_rmsd_between_poses_soluble>` for
    soluble proteins.

    Clustering is implemented by `SciPy <https://docs.scipy.org/doc/scipy-0.14.0/reference/cluster.hierarchy.html>`_.

    Parameters
    ----------
    dist_mat : np.array
        MxM distance matrix containing the RMSD between all pairs of
        poses. Where M is the number of poses.
    cutoff : float, optional
        The cutoff used to determine when to stop combining clusters.
        This happens when the shortest distance between all pairs of
        clusters is >cutoff. Default cutoff is 12.0 Angstroms.
    criterion : str, optional
        Default criterion is 'distance' which uses the cophenetic
        distance to compute distances between members from different
        clusters. The cophenetic distance is the distance between the 2
        clusters. See `SciPy documentation <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster>`_ for details.

    Returns
    -------
    np.array
        Returns an array of length M where each element in the array
        corresponds to cluster number that the m'th pose has been
        assigned to.
    """
    average_tree = sp.cluster.hierarchy.average(
        sp.spatial.distance.squareform(dist_mat))

    return sp.cluster.hierarchy.fcluster(average_tree, cutoff, criterion=criterion)


def calc_rot_static_coord(static_coord, mobile_com, static_com, rot_mat):
    """Compute coordinate of the static monomer after superimposition.

    After the mobile monomer is superimposed onto one of the monomers
    from the other pose, the static monomer moves with the mobile
    monomer. This function computes the coordinates of this static
    monomer.

    Parameters
    ----------
    static_coord : np.array
        Nx3 coordinate matrix of the static monomer with N atoms.
    mobile_com : np.array
        Vector pointing from the origin to the mobile monomer's center
        of mass (COM). COM of the mobile monomer in the pose before
        superimposition.
    static_com : np.array
        Vector pointing from the origin to the static monomer's center
        of mass (COM). COM of the static monomer in the pose before
        superimposition.
    rot_mat : np.array
        3x3 rotation matrix used to superimpose the mobile monomer to
        the monomer on the other pose.

    Returns
    -------
    np.array
        Nx3 coordinate matrix of the static monomer after superimposition
        of the mobile monomer.
    """
    return (static_coord - mobile_com).dot(rot_mat.T) + static_com


def cal_rmsd_between_poses_membrane(static_coord, mobile_coords):
    """Compute the RMSD/distance between the docking poses of HoTPs.

    Each pose contains 2 monomers, the static monomer whose coordinates
    are fixed and the mobile monomer that is translated and rotated
    during docking.
    
    There are 4 ways to compute the RMSD between any 2 monomers, each
    chosen from one of the 2 poses:

    1. Static monomers from both poses are superimposed and RMSD is computed for the mobile monomers. No superimposition is required for ZDOCK generated poses since the static monomers are fixed.
    2. Mobile monomers from both poses are superimposed and RMSD is computed for the static monomers.
    3. A mobile monomer from one pose is superimposed onto the static monomer in the second pose. The RMSD of monomers that were not superimposed are computed.
    4. The monomers used to compute the RMSD in (3) is superimposed and while the other two monomers are used to compute the RMSD.

    This function computes all 4 RMSDs and picks the smallest RMSD to
    represent the distance between these 2 poses.

    Parameters
    ----------
    static_coord : np.array
        Nx3 coordinate matrix of the static monomer with N atoms.
    mobile_coords : np.array
        MxNx3 coordinate matrices of the mobile monomers with N atoms
        from M poses.

    Returns
    -------
    np.array
        MxM matrix containing the RMSDs between the poses.
    """
    num_poses = mobile_coords.shape[0]
    num_atoms = mobile_coords.shape[1]

    static_com = static_coord.mean(axis=0)
    static_ori_coord = static_coord - static_com

    poses_mobile_com = mobile_coords.mean(axis=1)
    poses_mobile_ori_coords = mobile_coords - \
        poses_mobile_com[:, None, :]

    aln_mob_to_stat_static_coords = np.zeros(
        (num_poses, num_atoms, 3), dtype=np.float32)
    aln_mob_to_mob_static_coord = np.zeros((1, num_atoms, 3), dtype=np.float32)

    # Superimpose 2 poses at the static monomer and compute the MSD
    # between the 2 mobile monomers. Store them in this matrix:
    stat_stat_msd = np.zeros((num_poses, num_poses), dtype=np.float32)

    # Superimpose 2 poses at the mobile monomer and compute the MSD
    # between the 2 static monomers. Store them in this matrix:
    mob_mob_msd = np.zeros((num_poses, num_poses), dtype=np.float32)

    # Superimpose 2 poses at the mobile monomer and static monomer and
    # compute the MSD between the 2 remaining monomers. Store them in
    # this matrix:
    mob_stat_msd = np.zeros((num_poses, num_poses), dtype=np.float32)

    for row_idx in range(num_poses):
        # Each mobile pose vs the other mobile poses.
        # MSD matrix is symmetric therefore computing only the upper
        # half of the matrix.
        stat_stat_msd[row_idx, row_idx+1:] = calc_MSD(
            mobile_coords[row_idx], mobile_coords[row_idx+1:])

        # Superimpose the mobile monomer (pose 1) to the static monomer
        # (pose 2) and store the coordinates of the static monomer
        # (pose 1) after the superimposition.
        rot_mat_mobile_to_static = drsip_common.get_best_fit_rot_mat(
            poses_mobile_ori_coords[row_idx], static_ori_coord).astype('float32')
        aln_mob_to_stat_static_coords[row_idx, :] = calc_rot_static_coord(
            static_coord, poses_mobile_com[row_idx], static_com,
            rot_mat_mobile_to_static)

        for col_idx in range(row_idx+1, num_poses):
            # Superimpose the 2 poses at the mobile monomers and store
            # the static monomer after superimposition.
            rot_mat_mobile_to_mobile = drsip_common.get_best_fit_rot_mat(
                poses_mobile_ori_coords[row_idx],
                poses_mobile_ori_coords[col_idx]).astype('float32')
            aln_mob_to_mob_static_coord[0, :] = calc_rot_static_coord(
                static_coord, poses_mobile_com[row_idx],
                poses_mobile_com[col_idx], rot_mat_mobile_to_mobile)

            # Compute MSD between the 2 poses's static monomers for the
            # case where the mobile monomers are used for
            # superimposition. Symmetric matrix therefore we compute
            # only the upper half.
            mob_mob_msd[row_idx, col_idx] = calc_MSD(
                static_coord, aln_mob_to_mob_static_coord)

    # Compute the MSD between the mobile (pose 1) and static (pose 2)
    # monomer from the 2 poses after superimposition of the static
    # (pose 1) and mobile (pose 2) monomer and vice versa.
    # This matrix is not symmetric!
    mob_stat_msd[:] = calc_MSD_mob_static(
        aln_mob_to_stat_static_coords, mobile_coords)

    # Identify the minimum MSD between each pair of poses
    current_min_msd = min_dist(stat_stat_msd, mob_mob_msd)
    current_min_msd = min_dist(current_min_msd, np.triu(mob_stat_msd))
    current_min_msd = min_dist(current_min_msd, np.triu(mob_stat_msd.T))

    return np.sqrt(current_min_msd + current_min_msd.T)


def cal_rmsd_between_poses_soluble(static_coord, poses_mobile_coords):
    """Compute the RMSD between poses of soluble proteins.

    Parameters
    ----------
    static_coord : np.array
        Nx3 coordinate matrix of the static monomer with N atoms.
    mobile_coords : np.array
        MxNx3 coordinate matrices of the mobile monomers with N atoms
        from M poses.

    Returns
    -------
    np.array
        MxM matrix containing the RMSDs between the poses.
    """
    num_poses = poses_mobile_coords.shape[0]
    msd_between_conf = np.zeros((num_poses, num_poses), dtype=np.float32)

    for row_idx in xrange(num_poses):
        row_mobile_coord = poses_mobile_coords[row_idx]

        msd_between_conf[row_idx, row_idx +
                         1:] = calc_MSD(row_mobile_coord,
                                        poses_mobile_coords[row_idx+1:])

    return np.sqrt((msd_between_conf.T + msd_between_conf))
