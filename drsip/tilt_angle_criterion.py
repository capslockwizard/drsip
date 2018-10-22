"""
==================================================================
Tilt Angle Criterion (:mod:`drsip.tilt_angle_criterion`)
==================================================================

Module contains functions implementing the tilt angle criterion.

Functions
---------

.. autofunction:: tilt_angle_criterion
.. autofunction:: get_DSSP_assignment
.. autofunction:: determine_transmembrane_helices
.. autofunction:: get_cylindrical_axis
"""
import numpy as np
import pandas as pd
import Bio.PDB
import tempfile
import os
import re


def tilt_angle_criterion(cylindrical_axis, membrane_norm, cutoff=0.52359877559):
    """Tilt angle criterion

    Checks if the tilt angle of the monomer's cylindrical axis with the
    membrane normal. If the tilt angle is less than the cutoff then it
    passes the criterion otherwise it fails the criterion.

    Parameters
    ----------
    cylindrical_axis, membrane_norm : np.array
        3-dim vector of the monomer's cylindrical axis and axis that is
        orthogonal to the membrane bilayer.
    cutoff : float, optional
        Tilt angle cutoff in radians. Default cutoff: 0.52359877559 rad
        (30 degrees).

    Returns
    -------
    list
        Returns the tilt angle (float) and pass status (bool).
    """
    current_angle = np.arccos(cylindrical_axis.dot(membrane_norm))

    if current_angle > np.pi/2:
        current_angle = np.pi - current_angle

    pass_status = current_angle < cutoff

    return [current_angle, pass_status]


def get_DSSP_assignment(pdb_file, dssp_exec):
    """Obtain secondary structure assignment from DSSP.

    Run DSSP and obtain the secondary structure assignment of one chain
    in the PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file.
    dssp_exec : str
        Location of the DSSP executable

    Returns
    -------
    list
        Returns the DSSP results table (Dataframe) and the chain (str)
        that corresponds to the monomer submitted to DSSP. The DSSP
        results table contains the DSSP index, residue type,
        secondary structure type, solvent accessibility, backbone phi
        and psi angles.
    """
    idx = pd.IndexSlice

    # Use DSSP to assign secondary structure type to each residue.
    pdb_parser = Bio.PDB.PDBParser()
    structure = pdb_parser.get_structure('PDB', pdb_file)
    model = structure[0]
    dssp_results = Bio.PDB.DSSP(model, pdb_file, dssp=dssp_exec)

    first_chain = model.get_chains().next().get_id()

    chain_ID_resid_list = dssp_results.keys()
    dssp_results_table = pd.DataFrame(list(dssp_results))
    dssp_results_table = dssp_results_table[[0, 1, 2, 3, 4, 5]]
    dssp_results_table.columns = ['DSSP Index',
                                  'Residue Type', 'SS', 'Acc', 'Phi', 'Psi']
    dssp_results_table.index = pd.MultiIndex.from_tuples(
        [(chain_ID_resid_list[x][0], chain_ID_resid_list[x][1][1])
         for x in range(len(chain_ID_resid_list))])
    dssp_results_table = dssp_results_table.reindex()

    # Make sure all helices have negative phi and psi angles, otherwise
    # change them to '-'.
    helix_only_rows = (dssp_results_table['SS'] == 'H')
    helix_only_table = dssp_results_table[helix_only_rows]
    first_chain_helix_only = helix_only_table.loc[first_chain, :]
    pos_phi_or_psi_idx = first_chain_helix_only[(
        first_chain_helix_only['Phi'] >= 0).values +
        (first_chain_helix_only['Psi'] >= 0).values].index.tolist()
    dssp_results_table.loc[idx[first_chain, pos_phi_or_psi_idx], 'SS'] = '-'

    return [dssp_results_table, first_chain]


def determine_transmembrane_helices(sel, dssp_exec, min_res=15):
    """Determine the transmembrane helices

    The criterion for identifying transmembrane helices is: At least 15
    (min_res) consecutive residues are alpha-helices, allowing for at
    most 1 non-alpha helical residue in the middle of these residues.

    Parameters
    ----------
    sel : MDAnalysis.core.groups.AtomGroup
        MDAnalysis' atomgroup containing all the atoms of the monomer.
    dssp_exec : str
        Location of the DSSP executable.

    Returns
    -------
    list
        Returns the DSSP results table (Dataframe) and the chain (str)
        that corresponds to the monomer submitted to DSSP. The DSSP
        results table contains the DSSP index, residue type,
        secondary structure type, solvent accessibility, backbone phi
        and psi angles.
    """
    temp_pdb_file = tempfile.mkstemp(suffix='.pdb', prefix='DSSP', text=True)
    os.close(temp_pdb_file[0])
    sel.write(temp_pdb_file[1], file_format='pdb')
    dssp_results_table, first_chain = get_DSSP_assignment(
        temp_pdb_file[1], dssp_exec)
    os.remove(temp_pdb_file[1])

    # Given the secondary structure assignment, grab all non-overlapping helices
    ss_string = ''.join(dssp_results_table.loc[first_chain, 'SS'].tolist())
    ss_idx_resid = dssp_results_table.loc[first_chain, 'SS'].index.values

    ss_helix_element_pattern = re.compile('H+.{0,1}H+')
    helical_elements = ss_helix_element_pattern.finditer(ss_string)
    helical_elements_sel_str = []

    # Go through each helical element and make sure the length is at
    # least 15 residues long and that there are no missing residues in
    # the helix.
    #
    # Missing residues are identified by making sure that resids are
    # consecutive (no jumps).
    for helical_element in helical_elements:
        start_pos = helical_element.start()
        end_pos = helical_element.end()

        if (end_pos - start_pos) < min_res:
            continue

        helix_resids = ss_idx_resid[start_pos:end_pos]

        if np.where((helix_resids[1:] - helix_resids[:-1]) - 1)[0].size > 0:
            print 'Warning: Jump in resid in helix: %d-%d: %s' % (
                start_pos, end_pos, helix_resids)
            continue

        if helix_resids[-1] > helix_resids[0]:
            helical_elements_sel_str.append(
                'resid %d-%d' % (helix_resids[0], helix_resids[-1]))
        else:
            print 'Error: End resid is smaller than start resid: %d - %d' % (
                start_pos, end_pos)

    return [helical_elements_sel_str, dssp_results_table]


def get_cylindrical_axis(CA_sel, transmembrane_helix_sel_strs):
    """Get the cylindrical axis of a monomer

    The cylindrical axis of a monomer is the sum of the longitudinal
    axis (unit vector) of each transmembrane helix. The cylindrical
    axis is then normalized.

    Parameters
    ----------
    CA_sel : MDAnalysis.core.groups.AtomGroup
        MDAnalysis' atomgroup containing the monomer's C-alpha atoms.
    transmembrane_helix_sel_strs : list
        List of atom selection strings for each transmembrane helix.

    Returns
    -------
    np.array
        3-dim vector containing the cylindrical axis of the monomer.
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
        raise Exception('Selection is not a single chain')
