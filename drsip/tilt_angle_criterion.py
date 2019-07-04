"""
==================================================================
Tilt Angle Criterion (:mod:`drsip.tilt_angle_criterion`)
==================================================================

Module contains functions implementing the tilt angle criterion.

Functions
---------

.. autofunction:: tilt_angle_criterion
.. autofunction:: get_cylindrical_axis
"""
import numpy as np


def tilt_angle_criterion(cylindrical_axis, membrane_norm, cutoff=0.541052):
    """Tilt angle criterion

    The tilt angle of a monomer is the acute angle between its
    cylindrical axis and the membrane normal. If the tilt angle is less
    than the cutoff, it passes the criterion otherwise it fails the
    criterion.

    Parameters
    ----------
    cylindrical_axis, membrane_norm : np.array
        3-dim vector of the monomer's cylindrical axis and the membrane
        normal. The membrane normal is the axis orthogonal to the
        membrane bilayer.
    cutoff : float, optional
        Tilt angle cutoff in radians. Default cutoff is 0.541052 rad
        (31 degrees).

    Returns
    -------
    list
        Returns the tilt angle (float) and pass status (bool). The pass
        status is True if the monomer passes the filter, otherwise
        False.
    """
    current_angle = np.arccos(cylindrical_axis.dot(membrane_norm))

    if current_angle > np.pi/2:
        current_angle = np.pi - current_angle

    pass_status = current_angle < cutoff

    return [current_angle, pass_status]


def get_cylindrical_axis(CA_sel, transmembrane_helix_sel_strs):
    """Get the cylindrical axis of a monomer

    The cylindrical axis of a monomer is the average of the
    longitudinal axis (unit vector) of the transmembrane helices.

    Parameters
    ----------
    CA_sel : MDAnalysis.core.groups.AtomGroup
        MDAnalysis' atomgroup containing the monomer's C-alpha atoms.
    transmembrane_helix_sel_strs : list
        List of MDAnalysis style atom selection strings, selecting
        each transmembrane helix.

    Returns
    -------
    np.array
        3-dim vector (unit vector) containing the cylindrical axis of
        the monomer.
    """
    prot_segments = np.unique(CA_sel.segids)
    combined_transmembrane_helix_sel_str = '(' + \
        ' or '.join(transmembrane_helix_sel_strs) + ')'
    static_cylindrical_axis = np.zeros(3)

    first_static_helix_PA1 = None

    if len(prot_segments) == 1:
        transmembrane_helix_sel = CA_sel.select_atoms(
            combined_transmembrane_helix_sel_str)

        for transmembrane_helix_sel_str in transmembrane_helix_sel_strs:
            helix_sel = transmembrane_helix_sel.select_atoms(
                transmembrane_helix_sel_str)
            diff = helix_sel.positions - helix_sel.positions.mean(axis=0)

            # We do not scale the diff with 1/(N-1) because we're not using the
            # eigenvalues other than comparing their relative values (largest)
            static_helix_PA1 = np.linalg.eigh(diff.T.dot(diff))[1][:, 2]

            if first_static_helix_PA1 is None:
                first_static_helix_PA1 = static_helix_PA1

            elif static_helix_PA1.dot(first_static_helix_PA1) < 0:
                static_helix_PA1 = -static_helix_PA1

            static_cylindrical_axis += static_helix_PA1

        return static_cylindrical_axis/np.linalg.norm(static_cylindrical_axis)

    else:
        raise Exception('Residues from multiple chains were selected')
