"""
========================================================
DR-SIP package (:mod:`drsip`)
========================================================

Package contains filters and protocols for predicting homo-oligomeric
transmembrane proteins (HoTPs).

Functions
---------

.. autofunction:: DR_SIP
.. autofunction:: load_data

Classes
-------

.. autoclass:: DR_SIP_Base
    :members:

.. autoclass:: DR_SIP_Membrane
    :members:

.. autoclass:: DR_SIP_Soluble
    :members:

"""

from zdock_parser import ZDOCK
import pandas as pd
import numpy as np
import os
import string
from collections import Counter
import docking_eval
import MDAnalysis
import umsgpack
import zlib
import Cn_symm_criterion
import dist_restraints_filter
import struct_clustering
import tilt_angle_criterion
import save_load
import drsip_common as common


def load_data(filename):
    """
    Load saved DR_SIP data from file and return DR_SIP instance.

    Parameters
    ----------
    filename : str
        Full path to save file.

    Returns
    -------
    drsip_inst
        Instance of DR_SIP with data from save file.
    """
    with open(filename, 'rb') as input_file:
        temp_storage_dict = umsgpack.unpackb(
            zlib.decompress(input_file.read()))

    if temp_storage_dict['protocol'] == 'soluble':
        soluble_state = True

    elif temp_storage_dict['protocol'] == 'membrane':
        soluble_state = False

    else:
        raise ValueError('Invalid protocol type: {0}'.format(
            temp_storage_dict['protocol']))

    if temp_storage_dict['has_dist_rest_data']:
        save_load.load_StrIO(temp_storage_dict, 'dist_rest_file')
        dist_rest_file_input = temp_storage_dict.pop('dist_rest_file')
        save_load.load_pd_table(
            temp_storage_dict, 'dist_rest_conform_dist', pop=True)

    else:
        dist_rest_file_input = ''

    for file_var in ['static_pdb', 'mobile_pdb', 'zdock_output_file']:
        save_load.load_StrIO(temp_storage_dict, file_var)

    if not soluble_state:
        transmem_helix_sel_strs = temp_storage_dict.pop('transmem_helix_sel_strs')
    
    else:
        transmem_helix_sel_strs = []

    drsip_inst = DR_SIP(temp_storage_dict['static_pdb'], temp_storage_dict['mobile_pdb'],
                        temp_storage_dict['zdock_output_file'], dist_rest_file_input, transmem_helix_sel_strs, soluble_state)

    for pd_var in ['filter_data', 'filter_pass_results', 'filtered_poses_data', 'final_results']:
        save_load.load_pd_table(temp_storage_dict, pd_var, pop=True)

    if 'poses_CAPRI_class_val' in temp_storage_dict:
        save_load.load_pd_table(
            temp_storage_dict, 'poses_CAPRI_class', pop=True)

    if not soluble_state:
        temp_storage_dict['static_orient_unit_vec'] = np.array(
            temp_storage_dict['static_orient_unit_vec'], dtype=np.float32)
        temp_storage_dict['trans_vect'] = np.array(
            temp_storage_dict['trans_vect'], dtype=np.float32)
        temp_storage_dict['axis_of_rot'] = np.array(
            temp_storage_dict['axis_of_rot'], dtype=np.float32)
        temp_storage_dict['cylindrical_axis'] = np.array(
            temp_storage_dict['cylindrical_axis'], dtype=np.float32)

    drsip_inst.__dict__.update(temp_storage_dict)

    return drsip_inst


def DR_SIP(static_pdb, mobile_pdb, zdock_output_file, dist_rest_file='', transmem_helix_sel_strs=[], soluble=False):
    """
    Initialize the correct docking protocol.

    The default protocol is the membrane protocol. To select the
    soluble protocol, set soluble to True.

    Currently, only supports docking results from ZDOCK 3.0.2.

    Parameters
    ----------
    static_pdb, mobile_pdb : str
        Full path to the PDB files containing the static and mobile
        monomers.
    zdock_output_file : str
        Full path to the ZDOCK 3.0.2 generated output file.
    dist_rest_file : str, optional
        Full path to the distance restraints file, defaults to an empty
        string. When empty, the distance restraints filter is not used.
        Required when soluble is set to True.

        File is formatted such that each line corresponds to a distance
        restraint (tab delimited)::

            chain_id_1 res_id_1 chain_id_2 res_id_2 distance

    transmem_helix_sel_strs : lists of strings, optional
        Define which residues are the transmembrane helices. Each
        string should select for 1 transmembrane helix using the
        MDAnalysis' `selection language <https://www.mdanalysis.org/docs/documentation_pages/selections.html>`_. Required if soluble is False.
    soluble : bool, optional
        Defaults to False, which uses DR-SIP's membrane protocol. When
        set to True, uses the membrane protocol.
    """
    if soluble: # Soluble protocol

        if dist_rest_file != '':
            return DR_SIP_Soluble(static_pdb, mobile_pdb, zdock_output_file, dist_rest_file)

        else:
            raise IOError(
                'No distance restraints file')

    else: # Membrane protocol
        if len(transmem_helix_sel_strs) > 0: # Check if transmem_helix_sel_strs is provided
            return DR_SIP_Membrane(static_pdb, mobile_pdb, zdock_output_file, transmem_helix_sel_strs, dist_rest_file)
        
        else:
            raise ValueError("Transmembrane helix selection strings " \
            "(transmem_helix_sel_strs) are required to perform the " \
            "membrane protocol.")


class DR_SIP_Base(object):
    """
    DR_SIP base class.

    For subclassing only.

    Parameters
    ----------
    static_pdb, mobile_pdb : str
        Full path to the PDB files containing the static and mobile
        monomers.
    zdock_output_file : str
        Full path to the ZDOCK 3.0.2 generated output file.
    dist_rest_file : str, optional
        Full path to the distance restraints file, defaults to an empty
        string.

        File is formatted such that each line corresponds to a distance
        restraint (tab delimited)::

            chain_id_1 res_id_1 chain_id_2 res_id_2 distance

            chain_id_1 res_id_1 chain_id_2 res_id_2 distance
    """

    def __init__(self, static_pdb, mobile_pdb, zdock_output_file, dist_rest_file=''):
        self.static_pdb = static_pdb
        self.mobile_pdb = mobile_pdb
        self.zdock_output_file = zdock_output_file
        self.zdock_obj = ZDOCK(
            self.zdock_output_file, self.static_pdb, self.mobile_pdb)
        self.total_num_poses = self.zdock_obj.num_poses
        self.static_segids = np.unique(
            self.zdock_obj.static_uni.atoms.segids).tolist()
        self.mobile_segids = np.unique(
            self.zdock_obj.mobile_uni.atoms.segids).tolist()
        self.zdock_traj_uni = self.zdock_obj.get_MDAnalysis_Wrapper()

        if dist_rest_file != '':
            self.dist_rest_file = dist_rest_file
            self._parse_dist_res_input()
            self.has_dist_rest_data = True

        else:
            self.has_dist_rest_data = False

        self.run_status = False
        self.CAPRI_assign_status = False

    def _parse_dist_res_input(self):
        """Parses the distance restraints input file into a Pandas table."""
        self.dist_rest_table = pd.read_csv(
            self.dist_rest_file, delimiter='\t', header=None, index_col=None)
        self.dist_rest_table.columns = [
            'Segid 1', 'Resid 1', 'Segid 2', 'Resid 2', 'Dist']

        # Swap columns of rows to make sure that the segid 1 corresponds to
        # the static monomer
        swap_rows = np.where(
            ~self.dist_rest_table['Segid 1'].isin(self.static_segids))[0]
        self.dist_rest_table.iloc[swap_rows, 0:4] = \
            self.dist_rest_table.iloc[swap_rows][['Segid 2', 'Resid 2',
                'Segid 1', 'Resid 1']].values

        # Make sure each column's dtype is correct
        self.dist_rest_table = self.dist_rest_table.astype(
            dtype={'Segid 1': 'string', 'Resid 1': 'int32', 'Segid 2': 'string', 'Resid 2': 'int32', 'Dist': 'float32'})

        # Make sure the segids are valid segids
        segids_only = self.dist_rest_table.loc[:, [
            'Segid 1', 'Segid 2']]
        num_unknown_segids = np.where(~segids_only.isin(
            np.unique(self.static_segids + self.mobile_segids)))[0].shape[0]

        if num_unknown_segids > 0:
            raise Exception(
                'Invalid chain IDs in the distance restraints file: Must be from the static/mobile structure')

        # Check and see if there are pairs whose segids are identical
        num_identical_segids = np.where(
            segids_only.iloc[:, 0] == segids_only.iloc[:, 1])[0].shape[0]

        if num_identical_segids > 0:
            raise Exception(
                'Distance restraints file contains pairs of residues with the same chain ID')

        # Select of residues involved in the the distance restraints
        sel_str_format = '(segid {0} and resid {1})'.format

        static_sel_list = [sel_str_format(
            *resid) for resid in self.dist_rest_table[['Segid 1', 'Resid 1']].values]
        self.dist_res_static_sel_str = ' or '.join(
            static_sel_list) + ' and name CA'
        self.dist_res_static_sel = self.zdock_traj_uni.select_atoms(
            self.dist_res_static_sel_str)

        mobile_sel_list = [sel_str_format(
            *resid) for resid in self.dist_rest_table[['Segid 2', 'Resid 2']].values]
        self.dist_res_mobile_sel_str = ' or '.join(
            mobile_sel_list) + ' and name CA'
        self.dist_res_mobile_sel = self.zdock_traj_uni.select_atoms(
            self.dist_res_mobile_sel_str)

        seg_res_1_pairs_ordered = zip(
            self.dist_res_static_sel.segids, self.dist_res_static_sel.resids)
        seg_res_2_pairs_ordered = zip(
            self.dist_res_mobile_sel.segids, self.dist_res_mobile_sel.resids)

        # Indices for extracting the relevent residue pairs from the distance
        # matrix
        self.dist_mat_idx = [[seg_res_1_pairs_ordered.index(tuple(x[1].tolist())) for x in self.dist_rest_table.loc[:, ['Segid 1', 'Resid 1']].iterrows()], [
            seg_res_2_pairs_ordered.index(tuple(x[1].tolist())) for x in self.dist_rest_table.loc[:, ['Segid 2', 'Resid 2']].iterrows()]]

        # Distance restraints data
        self.dist_rest = self.dist_rest_table.loc[:, 'Dist'].values

        # Initialize table to store the distances of the residue pairs for
        # each pose
        col_name = ['{0}-{1}-{2}-{3}'.format(*self.dist_rest_table.iloc[idx, 0:4])
                    for idx in xrange(self.dist_rest_table.shape[0])]
        self.dist_rest_conform_dist = pd.DataFrame(np.zeros(
            (self.zdock_obj.num_poses, self.dist_rest_table.shape[0])), columns=col_name, index=range(1, self.zdock_obj.num_poses + 1))

    def write_results_file(self, filename):
        """
        Write the top results to a CSV file.

        Parameters
        ----------
        filename : str
            Full path to save the comma-seperated (CSV) file.
        """
        raise NotImplementedError

    def build_final_results_table(self, clusters, filtered_poses_data):
        """
        Returns the final results table.

        Table contains the final ranking of the docking poses.

        Needs to be implemented by the child class.

        Parameters
        ----------
        clusters : list of int
            List containing the cluster number that each pose is in.
        filtered_poses_data : pandas.Dataframe
            Table containing the filter data and cluster number of
            each pose.

        Returns
        -------
        results_table : pandas.Dataframe
            Final results table.
        """
        raise NotImplementedError

    def get_final_results_table(self):
        """Returns the final results table.

        Returns
        -------
        pandas.Dataframe
            Table of the final results.
        """
        if not self.run_status:
            self.run()

        return self.final_results

    def calc_CAPRI_class(self, ref_pdb_file, static_sel_str, mobile_sel_str):
        """
        Classify each docking pose to one of four CAPRI classes.

        Each docking pose is compared to the reference pose in
        ref_pdb_file and assigned to one of the following classes:
        'Incorrect', 'Acceptable', 'Medium' or 'High'.

        Uses the CAPRI criteria to classify the poses. See
        :py:class:`docking_eval.CAPRI_Criteria`.

        Needs to be implemented by the child class.

        Parameters
        ----------
        ref_pdb_file : str
            Full file path to PDB file containing the reference pose.
        static_sel_str, mobile_sel_str : str
            MDAnalysis style selection string to select the atoms in
            ref_pdb_file that corresponds to the static and mobile
            monomers, respectively.

        Returns
        -------
        CAPRI_class : pandas.Dataframe
            Table of the docking poses and their respective
            classifications.
        """
        raise NotImplementedError

    def get_CAPRI_class_table(self):
        """Returns the CAPRI classification table.

        Returns a table containing the CAPRI classification of each
        docking pose.

        Run the calc_CAPRI_class (:py:meth:`membrane <drsip.DR_SIP_Membrane.calc_CAPRI_class>` or :py:meth:`soluble <drsip.DR_SIP_Soluble.calc_CAPRI_class>`) method before
        calling this method.
        
        Returns
        -------
        pandas.Dataframe
            Table of the docking poses with their respective
            classifications, % Native Contacts and iRMSD.
        """
        if not self.CAPRI_assign_status:
            raise Exception('Run the calc_CAPRI_class method before running this method')

        return self.poses_CAPRI_class

    def run(self):
        """
        Executes the docking protocol.

        Needs to be implemented by the child class.
        """
        raise NotImplementedError


class DR_SIP_Membrane(DR_SIP_Base):
    """
    DR-SIP's membrane protocol.

    Protocol applies the C\ :sub:`n` symmetry, tilt angle and
    distance restraints (optional) filters on docking poses.

    Currently, only supports docking results from ZDOCK 3.0.2.

    Parameters
    ----------
    static_pdb, mobile_pdb : str
        Full path to the PDB files containing the static and mobile
        monomers.
    zdock_output_file : str
        Full path to the ZDOCK 3.0.2 generated output file.
    transmem_helix_sel_strs : lists of strings
        Define which residues are the transmembrane helices. Each
        string should select for 1 transmembrane helix using the
        MDAnalysis' `selection language <https://www.mdanalysis.org/docs/documentation_pages/selections.html>`_.
    dist_rest_file : str, optional
        Full path to the distance restraints file, defaults to an empty
        string.

        File is formatted such that each line corresponds to a distance
        restraint (tab delimited)::

            chain_1 resid_1 chain_2 resid_2 distance
    """
    def __init__(self, static_pdb, mobile_pdb, zdock_output_file, transmem_helix_sel_strs, dist_rest_file=''):
        super(DR_SIP_Membrane, self).__init__(static_pdb,
                                              mobile_pdb, zdock_output_file, dist_rest_file)

        if self.has_dist_rest_data:
            filter_data_num_cols = 6
            filter_pass_results_num_cols = 4
            filter_data_col_names = ['DR Spearman Corr', 'DR Pearson Corr',
                                     'Order of Symmetry', 'Cn Symm RMSD', 'Tilt Angle', 'Filter Pass Status']
            filter_pass_status_col_names = ['DR Pass Status', 'Cn Symm Pass Status',
                                            'Tilt Angle Pass Status', 'All Pass Status']

        else:
            filter_data_num_cols = 4
            filter_pass_results_num_cols = 3
            filter_data_col_names = [
                'Order of Symmetry', 'Cn Symm RMSD', 'Tilt Angle', 'Filter Pass Status']
            filter_pass_status_col_names = [
                'Cn Symm Pass Status', 'Tilt Angle Pass Status', 'All Pass Status']

        self.filter_data = pd.DataFrame(np.zeros((self.total_num_poses, filter_data_num_cols), dtype='float32'), index=range(
            1, self.total_num_poses + 1), columns=filter_data_col_names)
        self.filter_data.index.name = 'ZDOCK Rank'

        self.filter_pass_results = pd.DataFrame(np.zeros((self.total_num_poses, filter_pass_results_num_cols), dtype='bool'), index=range(
            1, self.total_num_poses + 1), columns=filter_pass_status_col_names)
        self.filter_pass_results.index.name = 'ZDOCK Rank'

        self.trans_vect = np.zeros(
            (self.total_num_poses, 3), dtype='float32')
        self.axis_of_rot = np.zeros(
            (self.total_num_poses, 3), dtype='float32')
        self.cylindrical_axis = np.zeros(3, dtype='float32')

        self.transmem_helix_sel_strs = transmem_helix_sel_strs

    def build_final_results_table(self, clusters, filtered_poses_data):
        """
        Returns the final results table.

        Table contains the representative poses of each cluster ranked
        based on the size of the clusters (in descending order). For each
        representative pose the Cluster ID, Order of Symmetry,
        C\ :sub:`n` symmetry RMSD, Tilt Angle, Pearson and Spearman
        correlations are included in the table.

        The representative pose of each cluster are chosen based on 2
        criteria:

        1. The consensus order of symmetry within a cluster. Only
           poses that have the consensus order of symmetry are chosen.
        2. The pose with the lowest C\ :sub:`n` symmetry RMSD.

        Use :py:meth:`get_final_results_table <drsip.DR_SIP_Soluble.get_final_results_table>` to return the table.

        Parameters
        ----------
        clusters : list of int
            List containing the cluster number that each pose is in.
        filtered_poses_data : pandas.Dataframe
            Table containing the filter data and cluster number of
            each pose.
        """
        cluster_counts = Counter(clusters)
        final_table = pd.DataFrame()

        for cluster_num, counts in cluster_counts.most_common():

            if counts < 3:
                break

            cluster_members = filtered_poses_data[filtered_poses_data['Cluster ID'] == cluster_num].copy(
            )
            consensus_order_of_symm = cluster_members['Order of Symmetry'].value_counts(
            ).index.values[0]

            consensus_cluster_members = cluster_members[cluster_members['Order of Symmetry'] == consensus_order_of_symm].copy(
            )
            consensus_cluster_members['C-n Symm Rank'] = np.argsort(
                consensus_cluster_members['Cn Symm RMSD'].values)
            final_table = final_table.append(
                consensus_cluster_members.sort_values('C-n Symm Rank').iloc[0].copy(), sort=True)

        if final_table.shape[0] != 0:
            final_table.drop('C-n Symm Rank', axis=1, inplace=True)
            final_table['DRSIP Rank'] = range(1, final_table.shape[0]+1)
            final_table.index.name = 'ZDOCK Rank'

        self.final_results = final_table.copy()

    def calc_CAPRI_class(self, ref_pdb_file, static_sel_str, mobile_sel_str, num_mers):
        """
        Classify each docking pose to one of four CAPRI classes.

        Each docking pose is compared to the reference pose in
        ref_pdb_file and assigned to one of the following classes:
        'Incorrect', 'Acceptable', 'Medium' or 'High'.

        Uses the modified CAPRI criteria for HoTPs to classify the
        poses. See :py:meth:`docking_eval.assign_Cn_symm_CAPRI_class <docking_eval.assign_Cn_symm_CAPRI_class>`.

        The modified CAPRI criteria includes an additional criterion
        where the predicted size ("n") of the pose must be equal the
        size of the reference complex (num_mers).

        Parameters
        ----------
        ref_pdb_file : str
            Full file path to PDB file containing the reference pose.
        static_sel_str, mobile_sel_str : str
            MDAnalysis style selection string to select the atoms in
            ref_pdb_file that corresponds to the static and mobile
            monomer, respectively.
        num_mers : int
            The size of the original C\ :sub:`n` symmetric HoTP complex.

        Returns
        -------
        CAPRI_class : pandas.Dataframe
            Table of the docking poses with their respective
            classifications, % Native Contacts and iRMSD.
        """
        
        if not self.run_status:
            self.run()

        universe = MDAnalysis.Universe(ref_pdb_file)
        self.zdock_traj_uni.trajectory[0]

        self.poses_CAPRI_class = docking_eval.assign_Cn_symm_CAPRI_class(
            universe.atoms, static_sel_str, mobile_sel_str, self.zdock_traj_uni.atoms, num_mers, self.filter_data['Order of Symmetry'].values)

        self.CAPRI_assign_status = True

    def run(self, Cn_rmsd_cutoff=2.0, dist_rest_cutoff=0.3, tilt_angle_cutoff=0.541052, cluster_RMSD_cutoff=12.0):
        """
        Executes the membrane protocol.

        Parameters
        ----------
        Cn_rmsd_cutoff : float, optional
            The C\ :sub:`n` symmetry RMSD filter cutoff value (in
            Angstrom). Default is 2 Angstrom.
        dist_rest_cutoff : float, optional
            The distance restraints filter cutoff value (correlation).
            Default correlation 0.3.
        tilt_angle_cutoff : float, optional
            The tilt angle filter cutoff value (radians). Default
            angle is 0.541052 rad (31 degrees).
        cluster_RMSD_cutoff : float, optional
            Cutoff used to cluster the docking poses. Default cutoff is
            12 Angstrom.
        """
        # Get the selection, coordinates and COM of the static monomer
        self.static_CA_sel = self.zdock_obj.static_uni.select_atoms(
            'segid ' + ' '.join(self.static_segids) + ' and name CA')
        static_CA_coord = self.static_CA_sel.positions.astype('float64')
        static_CA_com = static_CA_coord.mean(axis=0)
        static_CA_ori_coord = static_CA_coord - static_CA_com

        # Get the selection of the mobile monomer
        self.mobile_CA_sel = self.zdock_traj_uni.select_atoms(
            'segid ' + ' '.join(self.mobile_segids) + ' and name CA')

        # Temp arrays to store the filter results
        if self.has_dist_rest_data:
            static_dist_rest_CA_coord = self.dist_res_static_sel.positions
            tmp_filter_status = np.zeros(
                (self.total_num_poses, 4), dtype='bool')
            tmp_filter_results = np.zeros(
                (self.total_num_poses, 6), dtype='float32')

        else:
            tmp_filter_status = np.zeros(
                (self.total_num_poses, 3), dtype='bool')
            tmp_filter_results = np.zeros(
                (self.total_num_poses, 4), dtype='float32')

        self.cylindrical_axis = tilt_angle_criterion.get_cylindrical_axis(
            self.static_CA_sel, self.transmem_helix_sel_strs)

        self.static_orient_unit_vec = static_CA_ori_coord[0] / \
            np.linalg.norm(static_CA_ori_coord[0])

        # Temp arrays to store mobile coords that have passed the filters
        poses_mobile_coords = np.zeros(
            (self.total_num_poses, static_CA_coord.shape[0], 3), dtype='float32')
        poses_mobile_coords_idx = 0

        for pose_idx, ts in enumerate(self.zdock_traj_uni.trajectory):

            # Distance Restraints (DR) filter:
            if self.has_dist_rest_data:
                mobile_dist_rest_CA_coord = self.dist_res_mobile_sel.positions
                current_distances, spearman_correl, pearson_correl, dist_rest_pass_status = dist_restraints_filter.dist_restraints_filter(
                    static_dist_rest_CA_coord, mobile_dist_rest_CA_coord, self.dist_mat_idx, self.dist_rest, cutoff=dist_rest_cutoff)

                self.dist_rest_conform_dist.values[pose_idx] = current_distances
                tmp_filter_results[pose_idx, 0] = spearman_correl
                tmp_filter_results[pose_idx, 1] = pearson_correl
                tmp_filter_status[pose_idx, 0] = dist_rest_pass_status

            # Get the selection, coordinates and COM of the mobile monomer
            mobile_CA_coord = self.mobile_CA_sel.positions.astype(
                'float64')
            mobile_CA_com = mobile_CA_coord.mean(axis=0)
            mobile_CA_ori_coord = mobile_CA_coord - mobile_CA_com
            diff_vect = mobile_CA_com - static_CA_com

            # SIP filter:
            # Cn symmetry criterion
            axis_of_rot, static_trans_vect, num_mers, RMSD, C_n_pass_status = Cn_symm_criterion.Cn_symm_criterion(
                static_CA_ori_coord, mobile_CA_ori_coord, diff_vect, Cn_rmsd_cutoff)

            # Tilt angle criterion
            self.axis_of_rot[pose_idx] = axis_of_rot
            self.trans_vect[pose_idx] = static_trans_vect

            tilt_angle, tilt_angle_pass_status = tilt_angle_criterion.tilt_angle_criterion(
                self.cylindrical_axis, axis_of_rot, cutoff=tilt_angle_cutoff)

            # Store filter results
            if self.has_dist_rest_data:
                tmp_filter_results[pose_idx, 2] = num_mers
                tmp_filter_results[pose_idx, 3] = RMSD
                tmp_filter_status[pose_idx, 1] = C_n_pass_status

                tmp_filter_results[pose_idx, 4] = np.degrees(tilt_angle)
                tmp_filter_status[pose_idx, 2] = tilt_angle_pass_status

                tmp_filter_status[pose_idx, 3] = dist_rest_pass_status * \
                    C_n_pass_status * tilt_angle_pass_status
                tmp_filter_results[pose_idx,
                                   5] = tmp_filter_status[pose_idx, 3]

                if tmp_filter_status[pose_idx, 3]:
                    poses_mobile_coords[poses_mobile_coords_idx] = mobile_CA_coord
                    poses_mobile_coords_idx += 1

            else:
                tmp_filter_results[pose_idx, 0] = num_mers
                tmp_filter_results[pose_idx, 1] = RMSD
                tmp_filter_status[pose_idx, 0] = C_n_pass_status

                tmp_filter_results[pose_idx, 2] = np.degrees(tilt_angle)
                tmp_filter_status[pose_idx, 1] = tilt_angle_pass_status

                tmp_filter_status[pose_idx, 2] = C_n_pass_status * \
                    tilt_angle_pass_status
                tmp_filter_results[pose_idx,
                                   3] = tmp_filter_status[pose_idx, 2]

                if tmp_filter_status[pose_idx, 2]:
                    poses_mobile_coords[poses_mobile_coords_idx] = mobile_CA_coord
                    poses_mobile_coords_idx += 1

        self.filter_data.values[:] = tmp_filter_results
        self.filter_pass_results.values[:] = tmp_filter_status
        self.filtered_poses_data = self.filter_data[self.filter_pass_results['All Pass Status'].values].copy(
        ).drop('Filter Pass Status', axis=1)

        # Cluster remaining poses
        filtered_poses_dist_mat = struct_clustering.cal_rmsd_between_poses_membrane(
            static_CA_coord.astype('float32'), poses_mobile_coords[:poses_mobile_coords_idx])

        clusters = struct_clustering.cluster_poses(
            filtered_poses_dist_mat, cutoff=cluster_RMSD_cutoff)
        self.filtered_poses_data['Cluster ID'] = clusters

        # Select cluster representatives, then rank the clusters and
        # build the final results table
        self.build_final_results_table(
            clusters, self.filtered_poses_data)

        self.run_status = True

    def write_results_file(self, filename):
        """
        Write the top results to a CSV file.

        Parameters
        ----------
        filename : str
            Full path to save the comma-seperated (CSV) file.
        """
        if not self.run_status:
            self.run()

        common.makedir(filename)

        self.final_results.to_csv(filename)

    def save_data(self, filename):
        """Save current session/instance to file.

        Parameters
        ----------
        filename : str
            Full file path to save to
        """
        if not self.run_status:
            self.run()

        zdock_output_file_str = save_load.convert_StrIO_or_file_to_str(
            self.zdock_output_file)
        static_pdb_file_str = save_load.convert_StrIO_or_file_to_str(
            self.static_pdb)
        mobile_pdb_file_str = save_load.convert_StrIO_or_file_to_str(
            self.mobile_pdb)

        temp_storage_dict = {'has_dist_rest_data': self.has_dist_rest_data, 'static_pdb': static_pdb_file_str, 'mobile_pdb': mobile_pdb_file_str, 'zdock_output_file': zdock_output_file_str,
                             'static_orient_unit_vec': self.static_orient_unit_vec.tolist(), 'trans_vect': self.trans_vect.tolist(), 'axis_of_rot': self.axis_of_rot.tolist(),
                             'cylindrical_axis': self.cylindrical_axis.tolist(), 'protocol': 'membrane', 'transmem_helix_sel_strs': self.transmem_helix_sel_strs}

        save_load.save_pd_table(
            temp_storage_dict, self.filter_data, 'filter_data')
        save_load.save_pd_table(
            temp_storage_dict, self.filter_pass_results, 'filter_pass_results')
        save_load.save_pd_table(
            temp_storage_dict, self.filtered_poses_data, 'filtered_poses_data')
        save_load.save_pd_table(
            temp_storage_dict, self.final_results, 'final_results')

        if self.has_dist_rest_data:
            temp_storage_dict['dist_rest_file'] = save_load.convert_StrIO_or_file_to_str(
                self.dist_rest_file)
            save_load.save_pd_table(
                temp_storage_dict, self.dist_rest_conform_dist, 'dist_rest_conform_dist')

        if hasattr(self, 'poses_CAPRI_class'):
            save_load.save_pd_table(
                temp_storage_dict, self.poses_CAPRI_class, 'poses_CAPRI_class')

        common.makedir(filename)

        with open(filename, 'wb') as output_file:
            output_file.write(zlib.compress(
                umsgpack.packb(temp_storage_dict), 9))

    def write_pose(self, pose_num, filename):
        """Write a C\ :sub:`n` symmetric complex to a PDB file.

        Given the original rank of a ZDOCK pose, write out the closest
        ideal C\ :sub:`n` symmetric complex to a PDB file.

        Parameters
        ----------
        pose_num : int
            The pose's original ZDOCK rank.
        filename : str
            Filename to store the complex in PDB format.
        """
        pose_idx = pose_num - 1
        order_of_symm = np.int32(
            self.filter_data['Order of Symmetry'].values[pose_idx])
        rot_mat = Cn_symm_criterion.get_rot_mat_arbitrary_axis(
            2 * np.pi / order_of_symm, self.axis_of_rot[pose_idx])
        current_pos = self.zdock_obj.static_uni.atoms.positions - \
            self.zdock_obj.static_uni.atoms.center_of_mass() + \
            self.trans_vect[pose_idx]
        num_atoms = self.zdock_obj.static_uni.atoms.n_atoms

        complex_uni = MDAnalysis.Merge(
            *([self.zdock_obj.static_uni.atoms] * order_of_symm))

        current_select = complex_uni.select_atoms(
            'bynum 1:%d' % (num_atoms * 1))
        current_select.segments.segids = 'A'
        current_select.positions = current_pos

        for current_subunit_id in range(2, order_of_symm + 1):
            prev_subunit_id = current_subunit_id - 1
            current_pos = current_pos.dot(rot_mat.T)

            current_select = complex_uni.select_atoms('bynum %d:%d' % (
                num_atoms * prev_subunit_id + 1, num_atoms * current_subunit_id))
            current_select.positions = current_pos
            current_select.segments.segids = string.uppercase[prev_subunit_id]

        common.makedir(filename)

        complex_uni.atoms.write(filename)

    def write_topN_poses(self, folder, num_poses=20):
        """Write the top-N C\ :sub:`n` symmetric complexes to PDB files.

        The top-N complexes are from the DR-SIP final results table.

        Parameters
        ----------
        folder : str
            Folder to write the PDB files of the top-N complexes
            Files are written out as Pose_1.pdb, Pose_2.pdb, ...,
            Pose_N.pdb.
        num_poses : int, optional
            The number of top-N poses to write out. Default: 20
        """
        if not self.run_status:
            self.run()

        if self.final_results.shape[0] < num_poses:
            num_poses = self.final_results.shape[0]

        for idx in xrange(num_poses):
            pose_num = self.final_results.index.values[idx]

            self.write_pose(pose_num, os.path.join(
                folder, 'Pose_%d.pdb' % (idx + 1)))


class DR_SIP_Soluble(DR_SIP_Base):
    """
    DR-SIP's soluble protocol.

    Protocol applies the distance restraints filter on docking poses.

    Currently, only supports docking results from ZDOCK 3.0.2.

    Parameters
    ----------
    static_pdb, mobile_pdb : str
        Full path to the PDB files containing the static and mobile
        monomers.
    zdock_output_file : str
        Full path to the ZDOCK 3.0.2 generated output file.
    dist_rest_file : str
        Full path to the distance restraints file.

        File is formatted such that each line corresponds to a distance
        restraint (tab delimited)::

            chain_1 resid_1 chain_2 resid_2 distance
    """
    def __init__(self, static_pdb, mobile_pdb, zdock_output_file, dist_rest_file):
        super(DR_SIP_Soluble, self).__init__(static_pdb,
                                             mobile_pdb, zdock_output_file, dist_rest_file)

        self.filter_data = pd.DataFrame(np.zeros((self.total_num_poses, 3), dtype='float32'), index=range(
            1, self.total_num_poses + 1), columns=['DR Spearman Corr', 'DR Pearson Corr', 'Filter Pass Status'])
        self.filter_data.index.name = 'ZDOCK Rank'

        self.filter_pass_results = pd.DataFrame(np.zeros((self.total_num_poses), dtype='bool'), index=range(
            1, self.total_num_poses + 1), columns=['DR Pass Status'])
        self.filter_pass_results.index.name = 'ZDOCK Rank'

    def run(self, dist_rest_cutoff=0.3, cluster_RMSD_cutoff=12.0):
        """
        Executes the soluble protocol.

        Parameters
        ----------
        dist_rest_cutoff : float, optional
            The distance restraints filter cutoff value (correlation).
            Default correlation 0.3.
        cluster_RMSD_cutoff : float, optional
            Cutoff used to cluster the docking poses. Default cutoff is
            12 Angstrom.
        """

        # Get the selection, coordinates and COM of the static monomer
        self.static_CA_sel = self.zdock_obj.static_uni.select_atoms(
            'name CA')
        static_CA_coord = self.static_CA_sel.positions
        static_dist_rest_CA_coord = self.dist_res_static_sel.positions

        # Temp arrays to store the filter results
        tmp_filter_status = np.zeros(
            (self.total_num_poses, 1), dtype='bool')
        tmp_filter_results = np.zeros(
            (self.total_num_poses, 3), dtype='float32')

        # Temp arrays to store mobile coords that have passed the filters
        poses_mobile_coords = np.zeros((self.total_num_poses, self.zdock_obj.mobile_uni.select_atoms(
            'name CA').positions.shape[0], 3), dtype='float32')
        poses_mobile_coords_idx = 0

        self.mobile_CA_sel = self.zdock_traj_uni.select_atoms(
            'segid ' + ' '.join(self.mobile_segids) + ' and name CA')

        for pose_idx, ts in enumerate(self.zdock_traj_uni.trajectory):

            # Get the coordinates and COM of the mobile monomer
            mobile_CA_coord = self.mobile_CA_sel.positions.astype(
                'float64')
            mobile_dist_rest_CA_coord = self.dist_res_mobile_sel.positions

            # Distance Restraints (DR) filter:
            current_distances, spearman_correl, pearson_correl, dist_rest_pass_status = dist_restraints_filter.dist_restraints_filter(
                static_dist_rest_CA_coord, mobile_dist_rest_CA_coord, self.dist_mat_idx, self.dist_rest, cutoff=dist_rest_cutoff)

            # Store filter results
            self.dist_rest_conform_dist.values[pose_idx] = current_distances
            tmp_filter_results[pose_idx, 0] = spearman_correl
            tmp_filter_results[pose_idx, 1] = pearson_correl
            tmp_filter_results[pose_idx, 2] = dist_rest_pass_status
            tmp_filter_status[pose_idx] = dist_rest_pass_status

            if tmp_filter_status[pose_idx]:
                poses_mobile_coords[poses_mobile_coords_idx] = mobile_CA_coord
                poses_mobile_coords_idx += 1

        self.filter_data.values[:] = tmp_filter_results
        self.filter_pass_results.values[:] = tmp_filter_status
        self.filtered_poses_data = self.filter_data[self.filter_pass_results['DR Pass Status'].values].copy(
        ).drop('Filter Pass Status', axis=1)

        # Cluster remaining poses
        self.filtered_poses_dist_mat = struct_clustering.cal_rmsd_between_poses_soluble(
            static_CA_coord, poses_mobile_coords[:poses_mobile_coords_idx])

        clusters = struct_clustering.cluster_poses(
            self.filtered_poses_dist_mat, cutoff=cluster_RMSD_cutoff)
        self.filtered_poses_data['Cluster ID'] = clusters

        # Select cluster representatives, rank the clusters and build
        # the final results table
        self.build_final_results_table(
            clusters, self.filtered_poses_data)

        self.run_status = True

    def write_results_file(self, filename):
        """
        Write the top results to a CSV file.

        Parameters
        ----------
        filename : str
            Full path to save the comma-seperated (CSV) file.
        """
        if not self.run_status:
            self.run()

        common.makedir(filename)

        self.final_results.to_csv(filename)

    def build_final_results_table(self, clusters, filtered_poses_data):
        """
        Build the final results table.

        Table contains the representative poses of each cluster ranked
        based on the size of the clusters (in descending order). For each
        representative pose the Cluster ID, the Pearson and Spearman
        correlations are included in the table.

        The representative pose of each cluster has the highest Spearman
        correlation. If there are more than one pose with the same
        highest correlation, the one with the highest Pearson
        correlation is chosen.

        Use :py:meth:`get_final_results_table <drsip.DR_SIP_Soluble.get_final_results_table>` to return the table.

        Parameters
        ----------
        clusters : list of int
            List containing the cluster number that each pose is in.
        filtered_poses_data : pandas.Dataframe
            Table containing the filter data and cluster number of
            each pose.
        """
        cluster_counts = Counter(clusters)
        final_table = pd.DataFrame()

        for cluster_num, counts in cluster_counts.most_common():

            if counts < 3:
                break

            cluster_members = filtered_poses_data[filtered_poses_data['Cluster ID'] == cluster_num].copy(
            )
            cluster_members_max_Spearman = cluster_members[np.abs(
                cluster_members['DR Spearman Corr'] - cluster_members['DR Spearman Corr'].max()) < 10**(-5)]

            if cluster_members_max_Spearman.shape[0] == 1:
                final_table = final_table.append(
                    cluster_members_max_Spearman, sort=True)

            else:
                final_table = final_table.append(
                    cluster_members_max_Spearman.loc[cluster_members_max_Spearman['DR Pearson Corr'].idxmax()], sort=True)

        if final_table.shape[0] != 0:
            final_table['DRSIP Rank'] = range(1, final_table.shape[0]+1)
            final_table.index.name = 'ZDOCK Rank'

        self.final_results = final_table.copy()

    def calc_CAPRI_class(self, ref_pdb_file, static_sel_str, mobile_sel_str):
        """
        Classify each docking pose to one of four CAPRI classes.

        Each docking pose is compared to the reference pose in
        ref_pdb_file and assigned to one of the following classes:
        'Incorrect', 'Acceptable', 'Medium' or 'High'.

        Uses the CAPRI criteria to classify the poses. See
        :py:meth:`docking_eval.assign_soluble_CAPRI_class <docking_eval.assign_soluble_CAPRI_class>`.

        Parameters
        ----------
        ref_pdb_file : str
            Full file path to PDB file containing the reference pose.
        static_sel_str, mobile_sel_str : str
            MDAnalysis style selection string to select the atoms in
            ref_pdb_file that corresponds to the static and mobile
            monomer, respectively.

        Returns
        -------
        CAPRI_class : pandas.Dataframe
            Table of the docking poses with their respective
            classifications, % Native Contacts and iRMSD.
        """
        universe = MDAnalysis.Universe(ref_pdb_file)
        self.zdock_traj_uni.trajectory[0]

        self.poses_CAPRI_class = docking_eval.assign_soluble_CAPRI_class(
            universe.atoms, static_sel_str, mobile_sel_str, self.zdock_traj_uni.atoms)

        self.CAPRI_assign_status = True

        return self.poses_CAPRI_class

    def save_data(self, filename):
        """
        Save current session to file.

        Parameters
        ----------
        filename : str
            Full file path to save to
        """
        if not self.run_status:
            self.run()

        zdock_output_file_str = save_load.convert_StrIO_or_file_to_str(
            self.zdock_output_file)
        static_pdb_file_str = save_load.convert_StrIO_or_file_to_str(
            self.static_pdb)
        mobile_pdb_file_str = save_load.convert_StrIO_or_file_to_str(
            self.mobile_pdb)

        temp_storage_dict = {'has_dist_rest_data': self.has_dist_rest_data, 'static_pdb': static_pdb_file_str,
                             'mobile_pdb': mobile_pdb_file_str, 'zdock_output_file': zdock_output_file_str, 'protocol': 'soluble'}

        save_load.save_pd_table(
            temp_storage_dict, self.filter_data, 'filter_data')
        save_load.save_pd_table(
            temp_storage_dict, self.filter_pass_results, 'filter_pass_results')
        save_load.save_pd_table(
            temp_storage_dict, self.filtered_poses_data, 'filtered_poses_data')
        save_load.save_pd_table(
            temp_storage_dict, self.final_results, 'final_results')

        temp_storage_dict['dist_rest_file'] = save_load.convert_StrIO_or_file_to_str(
            self.dist_rest_file)
        save_load.save_pd_table(
            temp_storage_dict, self.dist_rest_conform_dist, 'dist_rest_conform_dist')

        if hasattr(self, 'poses_CAPRI_class'):
            save_load.save_pd_table(
                temp_storage_dict, self.poses_CAPRI_class, 'poses_CAPRI_class')

        common.makedir(filename)

        with open(filename, 'wb') as output_file:
            output_file.write(zlib.compress(
                umsgpack.packb(temp_storage_dict), 9))

    def write_pose(self, pose_num, filename):
        """
        Write a pose to a PDB file.

        Parameters
        ----------
        pose_num : int
            The pose's original ZDOCK rank.
        filename : str
            Filename of PDB file to store the docking pose.
        """
        self.zdock_obj.write_pose(pose_num, filename, mobile_only=False)

    def write_topN_poses(self, folder, num_poses=20):
        """
        Write the top-N docking poses to PDB.

        Parameters
        ----------
        folder : str
            Folder to write the PDB files of the top-N docking poses
            ranked by DR-SIP. Files are written out as Pose_1.pdb,
            Pose_2.pdb, ..., Pose_N.pdb.
        num_poses : int, optional
            The number of top poses to write out. Default: 20
        """
        if not self.run_status:
            self.run()

        if self.final_results.shape[0] < num_poses:
            num_poses = self.final_results.shape[0]

        for idx in xrange(num_poses):
            pose_num = self.final_results.index.values[idx]

            self.write_pose(pose_num, os.path.join(
                folder, 'Pose_%d.pdb' % (idx + 1)))
