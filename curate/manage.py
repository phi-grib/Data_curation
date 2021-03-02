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
import tarfile
import time

from curate.util import utils

from typing import Tuple, Union

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

    return True, 'curation repository updated\n'

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

#### API functions

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

def get_creation_date(endpoint_path: str) -> str:
    """
        Returns the creation date of the selected endpoint dir in the curation repository.

        :param endpoint_path: complete path to endpoint directory.

        :return creation_date: string with the creation date of the endpoint dir.
    """

    import datetime

    creation_date = datetime.datetime.fromtimestamp(os.stat(endpoint_path).st_mtime).strftime("%d-%b-%Y")
    
    return creation_date

def action_list(curation_dir: str) -> Tuple[bool, str]:
    """
        In no argument is provided lists all endpoints present at the repository 
        otherwyse lists all files for the endpoint provided as argument.

        :param curation_dir: path to the endpoint in curation repo
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
            creation_date = get_creation_date(xpath)
            sys.stderr.write("\n{} {}\n".format(x, creation_date))
            
        sys.stderr.write("\nRetrieved list of curation endpoints from {}\n".format(rdir))

        return True, "{} endpoints found\n".format(num_curs)

    else:
        # if a path name is provided, list files
        base_path = utils.curation_tree_path(curation_dir)
        num_files = 0
        sys.stderr.write('Files found in curation endpoint {}:\n'.format(curation_dir))
        for x in os.listdir(base_path):
            num_files += 1
            xpath = os.path.join(base_path,x)
            creation_date = get_creation_date(xpath)
            sys.stderr.write("\n{} {}\n".format(x, creation_date))

        return True, "Endpoint {} has {} published versions\n".format(curation_dir, num_files)

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

def action_dir() -> Tuple[bool,Union[str,list]]:
    """
        Returns a list of curation endpoints and files

        :return bool:
        :return str:
        :return results:
    """

    # get curation repo path

    cur_path = pathlib.Path(utils.curation_repository_path())
    if cur_path.is_dir() is False:
        return False,  'Curation repository path does not exist. Please run "flame -c config".\n'

    # get directories in curation repo path

    dirs = [x for x in cur_path.iterdir() if x.is_dir()]
    results = []

    for directory in dirs:
        dir_dict = {}
        # I convert directory, which is a PosixPath object into a string
        directory_string = str(directory)
        dir_dict[directory_string] = os.listdir(directory)
        results.append(dir_dict)

    return True, results

def action_kill(curation_endpoint: str) -> Tuple[bool,str]:
    """
        Removes the endpoint tree described by the argument.

        :param curation_endpoint: path to curation endpoint in the repo.

        :return bool:
        :return str:
    """

    if not curation_endpoint:
        return False, 'Empty endpoint name\n'

    ndir = utils.curation_tree_path(curation_endpoint)

    if not os.path.isdir(ndir):
        return False, "Model {} not found\n".format(curation_endpoint)

    try:
        shutil.rmtree(ndir, ignore_errors=True)
    except:
        return False, "Failed to remove model {}\n".format(curation_endpoint)

    sys.stderr.write("Model {} removed\n".format(curation_endpoint))
    
    return True, "Model {} removed\n".format(curation_endpoint)

def action_export(curation_endpoint: str) -> Tuple[bool,str]:
    """
        Exports the whole curation endpoint tree indicated in the argument as a single
        tarball file with the same name.

        :param curation_endpoint: path to curation endpoint in the repo.
    """

    if not curation_endpoint:
        return False,  'Empty endpoint name\n'

    current_path = os.getcwd()
    exportfile = os.path.join(current_path,curation_endpoint+'.tgz')

    base_path = utils.curation_tree_path(curation_endpoint)

    if not os.path.isdir(base_path):
        return False, 'Unable to export, endpoint directory not found\n'

    # change to curation repository to tar the file from there
    os.chdir(base_path)

    itemend = os.listdir()
    itemend.sort()

    with tarfile.open(exportfile, 'w:gz') as tar:
        for iversion in itemend:
            if not os.path.isdir(iversion):
                continue
            tar.add(iversion)

    # return to current directory
    os.chdir(current_path)
    sys.stderr.write("Endpoint {} exported as {}.tgz\n".format(curation_endpoint,curation_endpoint))

    return True, "Endpoint {} exported as {}.tgz\n".format(curation_endpoint,curation_endpoint)