"""
    Code that instantiate arguments comming from command line inputs

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 23/02/2021, 13:16 PM
"""

import os
import sys

import curate.manage as manage

from typing import Tuple, Optional

from curate.util import utils

def curation_cmd(commnad_dict: dict) -> Optional[bool]:
    """
        Instantiate curate objectt using commnad_dict from argument parser.

        :param commnad_dict:
                - data_input: input file name to be processed
                - molecule_identifier: column name containing the molecule ID. Usually CAS is used
                - endpoint: curation endpoint name
                - structure_column: column name containing the SMILES string
                - metadata: column names for metadata processing (only for API)
                - separator: file separator if input file is a csv or a tsv
                - remove_problematic: boolean indicating the option of removing problematic structures or not
                - outfile_type: output file type: xlsx, csv, tsv, sdf or json
    """
    
    import curate.dataset_curation as datacur
    
    # safety check if curation endpoint exists
    output_dir = utils.curation_tree_path(commnad_dict['endpoint'])
    if not os.path.isdir(output_dir):
        sys.stderr.write("Endpoint name not found in model repository.\n")
        return
    
    # check of metadata

    if commnad_dict['metadata']:
        metadata_ = commnad_dict['metadata']
        if (commnad_dict['molecule_identifier'] in metadata_) or (commnad_dict['structure_column'] in metadata_):
            sys.stderr.write("datacur curate : metadata can't contain the ID nor the SMILES column names.\n")
            return
    else:
        metadata_ = None
    
    # call of curation functions

    curating = datacur.DataCuration(data_input=commnad_dict['data_input'], 
                                    molecule_identifier=commnad_dict['molecule_identifier'],
                                    structure_column=commnad_dict['structure_column'],
                                    output_dir=output_dir,
                                    endpoint=commnad_dict['endpoint'],
                                    metadata=metadata_,
                                    separator=commnad_dict['separator'],
                                    remove_problematic=commnad_dict['remove_problematic'])

    curating.curate_data()
    curating.write_output_curation_data()

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
    elif arguments.action == 'export':
        success, results = manage.action_export(arguments.endpoint)
    elif arguments.action == 'download':
        success, results = manage.action_curation_results(arguments)
    else: 
        success = False
        results = "Specified manage action is not defined\n"

    return success, results