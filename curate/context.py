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
                - curation_type: chooses wheter to perform the classical chemical curation ('chem') or 
                a more specific one for htt files ('htt')
    """
    
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
    if commnad_dict['curation_type'] == 'chem':
        chemical_curation(commnad_dict=commnad_dict, output_dir=output_dir, metadata_=metadata_)
    elif commnad_dict['curation_type'] == 'htt':
        high_througput_curation(commnad_dict=commnad_dict, output_dir=output_dir, metadata_=metadata_)

def chemical_curation(commnad_dict: dict, output_dir: str, metadata_: str):
    """
        Performs the chemical curation of the file.
        Checks SMILES and fixes them if possible and adds a substance type based on the structure.

        :param command_dict:
        :param output_dir:
        :param metadata_:
    """

    import curate.dataset_curation as datacur

    curating = datacur.DataCuration(data_input=commnad_dict['data_input'], 
                                    molecule_identifier=commnad_dict['molecule_identifier'],
                                    structure_column=commnad_dict['structure_column'],
                                    output_dir=output_dir,
                                    endpoint=commnad_dict['endpoint'],
                                    metadata=metadata_,
                                    separator=commnad_dict['separator'],
                                    remove_problematic=commnad_dict['remove_problematic'],
                                    curation_type=commnad_dict['curation_type'])

    curating.curate_data()
    curating.write_output_curation_data()

def high_througput_curation(commnad_dict: dict, output_dir: str, metadata_: str):
    """
        Performs the curation of htt input files alongside with a chemical curation of the SMILES

        :param command_dict:
        :param output_dir:
    """

    import curate.htt.htt_curation as httcur
    
    htt_cur = httcur.htt_curation(data_input=commnad_dict['data_input'], 
                                    molecule_identifier=commnad_dict['molecule_identifier'],
                                    structure_column=commnad_dict['structure_column'],
                                    output_dir=output_dir,
                                    endpoint=commnad_dict['endpoint'],
                                    metadata=metadata_,
                                    separator=commnad_dict['separator'],
                                    remove_problematic=commnad_dict['remove_problematic'],
                                    curation_type=commnad_dict['curation_type'])


    htt_cur.chemical_curation()
    htt_cur.write_chem_output_and_x_matrix()

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