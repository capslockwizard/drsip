"""
================================================================
C\ :sub:`n` Symmetry Criterion (:mod:`drsip.Cn_symm_criterion`)
================================================================

Module contains functions implementing the C\ :sub:`n` symmetry criterion.

Functions
---------

.. autofunction:: Cn_symm_criterion
.. autofunction:: get_angle_of_rotation
.. autofunction:: get_axis_of_rot
.. autofunction:: get_rot_mat_arbitrary_axis
.. autofunction:: get_trans_vect
.. autofunction:: get_num_mers
"""
import numpy as np
import drsip_common


def Cn_symm_criterion(static_ori_coord, mobile_ori_coord, diff_vect, cutoff=2.0):
    """
    C\ :sub:`n` symmetry criterion.

    Checks if the current docking pose is close to its closest ideal C\ :sub:`n`
    symmetric pose.

    Parameters
    ----------
    static_ori_coord, mobile_ori_coord : np.array
        Nx3 coordinate array of the static and mobile monomers,
        respectively, from the docking pose. Where N is the number of
        atoms.
    diff_vect : np.array
        Vector from the center of mass (COM) of the static monomer to
        the COM of the mobile monomer.
    cutoff : float, optional
        The C\ :sub:`n` symmetric RMSD cutoff. If the docking pose's
        C\ :sub:`n` symmetric RMSD > cuttoff, it fails the criterion.
        Default cutoff is 2.0 Angstrom.

    Returns
    -------
    list
        Returns the axis of rotation (np.array), static monomer's
        translation vector (np.array), number of oligomers (int), C\ :sub:`n`
        symmetric RMSD (float) and RMSD_pass_status (bool). The
        RMSD_pass_status is True when the docking pose passes the
        criterion, otherwise False.
    """
    rot_mat = drsip_common.get_best_fit_rot_mat(
        static_ori_coord, mobile_ori_coord)
    axis_of_rot = get_axis_of_rot(rot_mat)
    rot_angle = get_angle_of_rotation(rot_mat, axis_of_rot)

    num_mers = get_num_mers(rot_angle)
    ideal_rot_angle = np.sign(rot_angle) * 2 * np.pi / num_mers

    if num_mers == 2:
        diff_vect_ortho = diff_vect - \
            diff_vect.dot(axis_of_rot) * axis_of_rot
        static_trans_vect = -diff_vect_ortho / 2

    else:
        static_trans_vect = get_trans_vect(
            axis_of_rot, ideal_rot_angle, diff_vect)

    ideal_rot_mat = get_rot_mat_arbitrary_axis(
        ideal_rot_angle, axis_of_rot)

    RMSD_pass_status = True

    pred_mobile_coord = (static_ori_coord +
                         static_trans_vect).dot(ideal_rot_mat.T)
    RMSD = Cn_symm_RMSD(pred_mobile_coord,
                        mobile_ori_coord, static_trans_vect, diff_vect)

    if RMSD > cutoff:
        RMSD_pass_status = False

    return [axis_of_rot, static_trans_vect, num_mers, RMSD, RMSD_pass_status]


def Cn_symm_RMSD(pred_mobile_coord, mobile_ori_coord, static_trans_vect, diff_vect):
    """
    Calculate the C\ :sub:`n` symmetry RMSD.

    Parameters
    ----------
    pred_mobile_coord : np.array
        Nx3 coordinate array of the ideal C\ :sub:`n` symmetric mobile
        monomer. Where N is the number of atoms.
    mobile_ori_coord : np.array
        Nx3 coordinate array of the mobile monomer from the docking
        pose. Where N is the number of atoms.
    static_trans_vect : np.array
        The vector pointing from the center of mass (COM) of the
        complex to the COM of the static monomer. The COM of the
        complex is set at the origin.
    diff_vect : np.array
        Vector from the center of mass (COM) of the static monomer to
        the COM of the mobile monomer.

    Returns
    -------
    np.array
        Returns the axis of rotation.
    """
    return np.sqrt(np.sum((pred_mobile_coord - mobile_ori_coord -
                   static_trans_vect - diff_vect)**2) /
                   pred_mobile_coord.shape[0])


def get_axis_of_rot(rot_mat):
    """
    Extract the axis of rotation from the rotation matrix.

    Parameters
    ----------
    rot_mat : np.array
        3x3 rotation matrix.

    Returns
    -------
    np.array
        Returns the axis of rotation.
    """
    rot_mat_eigval, rot_mat_eigvect = np.linalg.eig(rot_mat)
    eigenvec_idx = np.argmin(np.abs(np.real(rot_mat_eigval) - 1))

    return np.float64(np.real(rot_mat_eigvect[:, eigenvec_idx]))


def get_angle_of_rotation(rot_mat, axis_of_rot):
    """
    Extract the angle of rotation from the rotation matrix.

    Parameters
    ----------
    rot_mat : np.array
        3x3 rotation matrix.
    axis_of_rot : np.array
        Axis of rotation extracted from the rotation matrix. See
        :py:meth:`get_axis_of_rot <drsip.Cn_symm_criterion.get_axis_of_rot>`.

    Returns
    -------
    float
        Returns the angle of rotation in radians.
    """
    skew_symm_mat = rot_mat - rot_mat.T

    sin_angle_rot = (skew_symm_mat[2, 1] + skew_symm_mat[0, 2] +
                     skew_symm_mat[1, 0]) / (2 * np.sum(axis_of_rot))
    cos_angle_rot = (np.trace(rot_mat) - 1) / 2

    return np.arctan2(sin_angle_rot, cos_angle_rot)


def get_num_mers(rot_angle):
    """
    Predict the size of the complex from the rotation angle.

    The size ("n") is the number of monomers in the C\ :sub:`n` symmetric
    complex.

    Parameters
    ----------
    rot_angle : float
        Rotation angle in radians. See :py:meth:`get_angle_of_rotation <drsip.Cn_symm_criterion.get_angle_of_rotation>`.

    Returns
    -------
    int
        Number of monomers ("n").
    """
    return np.abs(np.rint((np.pi * 2) / rot_angle))


def get_trans_vect(axis_of_rot, rot_angle, diff_vect):
    """
    Returns the translation vector.

    The translation vector is the vector from the complex's center of
    mass (COM) to the COM of the static monomer. The COM of the complex
    is set to the origin.

    Parameters
    ----------
    axis_of_rot : np.array
        Axis of rotation. See :py:meth:`get_axis_of_rot <drsip.Cn_symm_criterion.get_axis_of_rot>`.
    rot_angle : float
        Rotation angle in radians. See :py:meth:`get_angle_of_rotation <drsip.Cn_symm_criterion.get_angle_of_rotation>`.
    diff_vect : np.array
        Vector from the center of mass (COM) of the static monomer to
        the COM of the mobile monomer.

    Returns
    -------
    trans_vect : np.array
        Translation vector.
    """
    diff_vect_ortho = diff_vect - diff_vect.dot(axis_of_rot) * axis_of_rot
    dist_between_com = np.linalg.norm(diff_vect_ortho)
    dist_between_com_n_axis_rotation = np.abs(
        dist_between_com / (2 * np.sin(rot_angle / 2)))
    diff_vect_to_trans_vect_rot_mat = get_rot_mat_arbitrary_axis(
        -np.sign(rot_angle) * (np.abs(rot_angle) + np.pi) / 2, axis_of_rot)

    trans_vect = diff_vect_to_trans_vect_rot_mat.dot(diff_vect_ortho)
    trans_vect /= np.linalg.norm(trans_vect)
    trans_vect *= dist_between_com_n_axis_rotation

    return trans_vect


def get_rot_mat_arbitrary_axis(rot_angle, axis_of_rot):
    """
    Returns the rotation matrix for the rotation by the given rotation
    angle about the axis of rotation.

    Parameters
    ----------
    rot_angle : float
        Rotation angle in radians.
    axis_of_rot : np.array
        Axis of rotation.

    Returns
    -------
    rot_mat : np.array
        3x3 rotation matrix.
    """
    axis_of_rot = np.array(axis_of_rot)
    x, y, z = axis_of_rot / np.linalg.norm(axis_of_rot)

    rot_mat = np.zeros((3, 3))

    sin_theta = np.sin(rot_angle)
    cos_theta = np.cos(rot_angle)
    minus_cos_theta = 1 - cos_theta

    rot_mat[0, 0] = cos_theta + x**2 * minus_cos_theta
    rot_mat[0, 1] = x * y * minus_cos_theta - z * sin_theta
    rot_mat[0, 2] = x * z * minus_cos_theta + y * sin_theta
    rot_mat[1, 0] = x * y * minus_cos_theta + z * sin_theta
    rot_mat[1, 1] = cos_theta + y**2 * minus_cos_theta
    rot_mat[1, 2] = y * z * minus_cos_theta - x * sin_theta
    rot_mat[2, 0] = x * z * minus_cos_theta - y * sin_theta
    rot_mat[2, 1] = y * z * minus_cos_theta + x * sin_theta
    rot_mat[2, 2] = cos_theta + z**2 * minus_cos_theta

    return rot_mat
