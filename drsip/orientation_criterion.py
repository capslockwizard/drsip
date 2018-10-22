"""
===============================================================
Orientation Criterion (:mod:`drsip.orientation_criterion`)
===============================================================

Module contains the function implementing the orientation criterion.

Functions
---------

.. autofunction:: orientation_criterion
"""


def orientation_criterion(mobile_orient_unit_vec, static_orient_unit_vec, parallel=True):
    """Orientation criterion.

    The 2 monomers in the docking pose can be parallel or
    anti-parallel. The orientation of each monomer is defined as the
    vector from the center of mass of monomer to the N-terminus
    residue's C-alpha atom. The direction is parallel when the dot
    product of the 2 vectors from the 2 monomers is positive (<90
    degrees) otherwise they are anti-parallel.

    Determines if the orientation of the pose passes the criterion
    (parallel or anti-parallel).

    Parameters
    ----------
    mobile_orient_unit_vec, static_orient_unit_vec : np.array
        3-dim vector pointing from the center of mass the monomer to
        the N-term residue's C-alpha atom.
    parallel : bool, optional
        Check if the pose is parallel (default: True) or anti-parallel
        (False).

    Returns
    -------
    bool
        True if the pose passes the criterion and False otherwise.
    """

    if parallel:
        return (mobile_orient_unit_vec.dot(static_orient_unit_vec)) > 0
    else:
        return (mobile_orient_unit_vec.dot(static_orient_unit_vec)) <= 0
