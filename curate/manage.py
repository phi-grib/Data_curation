"""
    Managing functions for Data curation CLI.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/02/2021, 17:49 PM
"""

import os
import pandas as pd
import pathlib
import shutil
import sys
import time

from curate.util import utils

from typing import Tuple

def get_metadata(data: pd.DataFrame, structure_colname: str) -> list:
    """
        Checks the columns of the input data and returns all that are 
        not the SMILES. Everything but the structure is considered metadata.
        
        :param data: dataframe to be used
        :param structure_colname: name of the column containing the SMILES

        :return metadata: returns a list with the metadata, which are the 
        column names of the input data without the SMILES
    """

    metadata = data.loc[:, data.columns != structure_colname].columns

    return metadata

def convert_to_json(data: pd.DataFrame, filename: str = None) -> str:
    """
        Converts into JSON the input dataframe.
        If filename is provided, json is written as a file.

        :param data: dataframe to be converted.
        :param filename: filename where to dump the JSON.

        :return json_data: JSON string from pandas dataframe
    """

    json_data = data.to_json(path_or_buf = filename, orient = 'index')

    return json_data

def set_curation_repository(path: str = None):
    """
        Set the path to the curation repository

        :param path: string indicating the path of the curation repository
    """

    utils.set_curation_repository(path)

    sys.stderr.write('Model repository updated to {}'.format(path))

    return True, 'curation repository updated'

#### Functions to represent simple statistics after the curation

def get_number_of_processed_vs_unprocessed(smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
        This function returns a dataframe with the number of total SMILES, the ones
        that have been processed by the code and the ones that haven't.

        Is it expected to receive an input from the API using the output of the previous
        function get_list_of_filtered_SMILES_and_substance_type()

        :param smiles_dataframe:

        :return smiles_stats_df:
    """

    unprocessed_smiles = smiles_dataframe[smiles_dataframe['structure_curated'].isna()]
    processed_smiles = smiles_dataframe[~smiles_dataframe['structure_curated'].isna()]

    smiles_stats_dict = {'Total SMILES':len(smiles_dataframe.index), 
                         'Processed SMILES':len(processed_smiles.index), 
                         'Unable to process':len(unprocessed_smiles.index)}
    
    smiles_stats_df = pd.DataFrame(data=smiles_stats_dict)

    return smiles_stats_df

def get_total_of_smiles_per_type_of_substance(smiles_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
        This function returns the amount of substance types found of each kind
        after curating the SMILES.

        :param smiles_dataframe:

        :return subs_count:
    """

    subs_count = smiles_dataframe.groupby('substance_type_name')['substance_type_name'].count()

    return subs_count

def action_new(curation_path: str) -> Tuple[bool, str]:
    """
        Create a new curation endpoint tree, using the given name.
        
        :param curation_path: curation endpoint in curation repository where output will be saved

        :return bool: True when evertyhing has workded, otherwise False.
        :return str: strings that would be the equivalent to the standard error.
    """

    if not curation_path:
        return False, 'empty endpoint curation label\n'

    # importlib does not allow using 'test' and issues a misterious error when we
    # try to use this name. This is a simple workaround to prevent creating paths 
    # with this name 
    if curation_path == 'test':
        return False, 'the name "test" is disallowed, please use any other name\n'

    # curation endpoint directory
    ndir = pathlib.Path(utils.curation_tree_path(curation_path))

    # check if there is already a tree for this endpoint
    if ndir.exists():
        return False, "Endpoint {} already exists\n".format(curation_path)

    try:
        ndir.mkdir(parents=True)
        sys.stderr.write("{} created\n".format(ndir))
    except:
        return False, "Unable to create path for {} endpoint\n".format(curation_path)

    sys.stderr.write("New endpoint {} created\n".format(curation_path))
    
    return True, "new endpoint {} created\n".format(curation_path)

def action_list(curation_dir: str) -> Tuple[bool, str]:
    """
        In no argument is provided lists all models present at the repository 
        otherwyse lists all versions for the model provided as argument
    """

    # if no name is provided, just list the different curation dirs
    if not curation_dir:
        rdir = utils.curation_repository_path()
        if os.path.isdir(rdir) is False:
            return False, 'the curation repository path does not exist. Please run "datacur -c config".\n'

        num_curs = 0
        sys.stderr.write('Curation endpoints found in repository:\n')
        for x in os.listdir(rdir):
            xpath = os.path.join(rdir,x) 
            # discard if the item is not a directory
            if not os.path.isdir(xpath):
                continue
            num_curs += 1
            sys.stderr.write("\t"+x)
        sys.stderr.write("\nRetrieved list of curation endpoints from {}\n".format(rdir))
        return True, "{} endpoints found\n".format(num_curs)

def action_remove(curation_endpoint: str) -> Tuple[bool, str]:
    """
        Remove the curation endpoint directory indicated as 
        argument

        :param curation_endpoint: curation endpoint to be removed

    """

    if not curation_endpoint:
        return False, 'Empty curation endpoint\n'

    rdir = utils.curation_tree_path(curation_endpoint)
    if not os.path.isdir(rdir):
        return False, '{curation_endpoint} not found\n'.format(curation_endpoint)

    shutil.rmtree(rdir, ignore_errors=True)
    sys.stderr.write("Curation endpoint dir {} has been removed\n".format(curation_endpoint))

    return True, "Curation endpoint dir {} has been removed\n".format(curation_endpoint)
