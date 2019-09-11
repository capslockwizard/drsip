from zdock_parser import ZDOCK
import numpy as np
import save_load


class ZDOCK_Reader(object):

    def __init__(self, static_pdb, mobile_pdb, zdock_output_file):
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
        self.static_sel = self.zdock_obj.static_uni.atoms
        self.mobile_sel = self.zdock_obj.mobile_uni.atoms
        self.poses_uni = self.zdock_obj.get_MDAnalysis_Wrapper()

    def save_data(self, storage_dict):
        """Save current session/instance to file.

        Parameters
        ----------
        filename : str
            Full file path to save to
        """

        zdock_output_file_str = save_load.convert_StrIO_or_file_to_str(
            self.zdock_output_file)
        static_pdb_file_str = save_load.convert_StrIO_or_file_to_str(
            self.static_pdb)
        mobile_pdb_file_str = save_load.convert_StrIO_or_file_to_str(
            self.mobile_pdb)

        temp_storage_dict = {'static_pdb': static_pdb_file_str,
                             'mobile_pdb': mobile_pdb_file_str,
                             'zdock_output_file': zdock_output_file_str,
                             'reader': 'ZDOCK',
                             }

        storage_dict.update(temp_storage_dict)

        return storage_dict

    @staticmethod
    def load_data(storage_dict):

        for file_var in ['static_pdb', 'mobile_pdb', 'zdock_output_file']:
            save_load.load_StrIO(storage_dict, file_var)

        return ZDOCK_Reader(storage_dict['static_pdb'],
                            storage_dict['mobile_pdb'],
                            storage_dict['zdock_output_file']
                            )


class MDA_Reader(object):

    def __init__(self, static_sel, mobile_sel, poses_uni):
        self.static_sel = static_sel
        self.mobile_sel = mobile_sel
        self.poses_uni = poses_uni
        self.total_num_poses = self.poses_uni.trajectory.n_frames
        self.static_segids = np.unique(
            self.static_sel.segids).tolist()
        self.mobile_segids = np.unique(
            self.mobile_sel.segids).tolist()

    def save_data(self, storage_dict):
        """Save current session/instance to file.

        Parameters
        ----------
        filename : str
            Full file path to save to
        """

        storage_dict['reader'] = 'MDA'

        return storage_dict

    @staticmethod
    def load_data(storage_dict):
        return None


readers_dict = {'ZDOCK': ZDOCK_Reader, 'MDA': MDA_Reader}
