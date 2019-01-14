"""
=====================================================
Saving and Loading (:mod:`drsip.save_load`)
=====================================================

Module contains helper functions to save and load DR-SIP data files.

Functions
---------

.. autofunction:: save_pd_table
.. autofunction:: load_pd_table
.. autofunction:: load_StrIO
.. autofunction:: convert_StrIO_or_file_to_str
"""

import umsgpack
import zlib
import pandas as pd
import drsip_common


def save_pd_table(storage, table, var_name):
    """Split and store the values and indices 
    
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

    storage[var_name + '_val'] = table.values.tolist()
    storage[var_name + '_idx'] = table.index.tolist()

    if table.index.name is None:
        storage[var_name + '_col'] = table.reset_index().columns.tolist()

    else:
        storage[var_name + '_col'] = table.columns.tolist()


def load_pd_table(storage, var_name, pop=False):
    """Reconstruct Pandas table from dictionary key-value pairs and save to storage

    Args:
        storage (dict): Dictionary to load the data from
        var_name (str): Original variable name
        pop (bool, optional): If True, remove key-value pairs used to reconstruct the table from storage
    """

    if pop:
        value = storage.pop(var_name + '_val')
        index = storage.pop(var_name + '_idx')
        columns = storage.pop(var_name + '_col')

    else:
        value = storage[var_name + '_val']
        index = storage[var_name + '_idx']
        columns = storage[var_name + '_col']

    if len(value) != 0 and len(value[0]) != len(columns):
        storage[var_name] = pd.DataFrame(
            value, index=index, columns=columns[1:])
        storage[var_name].index.name = columns[0]

    else:
        storage[var_name] = pd.DataFrame(value, index=index, columns=columns)


def load_StrIO(storage, var_name):
    """Convert strings in storage to StringIO and replace original string in storage

    Args:
        storage (dict): Dictionary to load and save the string/StringIO data
        var_name (str): Variable name
    """

    file_str = storage.pop(var_name)

    storage[var_name] = drsip_common.convert_str_to_StrIO(file_str)


def convert_StrIO_or_file_to_str(input_obj):

    if isinstance(input_obj, str):
        input_str = drsip_common.convert_file_to_str(input_obj)

    else:
        input_str = drsip_common.convert_StrIO_to_str(input_obj)

    return input_str
