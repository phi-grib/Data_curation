"""
    Code that instantiate arguments comming from command line inputs

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 23/02/2021, 13:16 PM
"""
import os
import sys

import curate.manage as manage

from typing import Tuple

from curate.util import utils

def curation_cmd(**kwargs: dict):
    """
        Instantiate curate objectt using kwargs from argument parser.

        :param **kwargs:
                - data_input: input file name to be processed
                - molecule_identifier: column name containing the molecule ID. Usually CAS is used.
                - structure_column: column name containing the SMILES string
                - separator: file separator ir input file is a csv or a tsv.
                - remove_problematic: boolean indicating the option of removing problematic structures or not
                - outfile_name: output file name
                - outfile_type: output file type: xlsx, csv, tsv or sdf.
    """
    
    import curate.dataset_curation as datacur
    
    # safety check if curation endpoint exists
    output_dir = utils.curation_tree_path(kwargs['endpoint'])
    if not os.path.isdir(output_dir):
        sys.stderr.write("Endpoint name not found in model repository.\n")
        return

    curating = datacur.DataCuration(data_input=kwargs['data_input'], 
                                    molecule_identifier=kwargs['molecule_identifier'],
                                    structure_column=kwargs['structure_column'],
                                    output_dir=output_dir,
                                    separator=kwargs['separator'])

    curating.curate_data(remove_problematic=kwargs['remove_problematic'])

    curating.get_output_file(outfile_name=kwargs['outfile_name'],
                             outfile_type=kwargs['outfile_type'])

def manage_cmd(arguments: dict) -> Tuple[bool, str]:
    """
        Calls diverse maintenance commands

        :param arguments: arguments from argparser

        :return bool:
        :return str:
    """

    if arguments.action == 'new':
        success, results = manage.action_new(arguments.endpoint)
    elif arguments.action == 'list':
        success, results = manage.action_list(arguments.endpoint)
    elif arguments.action == 'remove':
        success, results = manage.action_remove(arguments.endpoint)
    else: 
        success = False
        results = "Specified manage action is not defined"

    return success, results

def output_file_to_proper_dir(outfile_name: str) -> str:
    """
        Join curation repo path to output file name so it can be saved 
        in its proper directory.

        :param outfile_name: output file name

        :return outfile_path: output file name with the directory path
    """
    
    curation_dir_path = utils.curation_repository_path()
    outfile_path = '/'.join([curation_dir_path,outfile_name])

    return outfile_path