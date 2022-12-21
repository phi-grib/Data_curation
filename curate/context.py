"""
    Code that instantiate arguments comming from command line inputs

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 23/02/2021, 13:16 PM
"""

import os
import sys

import curate.metacuration.metacuration as metacur
import curate.manage as manage

from typing import Tuple, Optional

from curate.util import utils, get_logger
LOG = get_logger(__name__)

def curation_cmd(command_dict: dict) -> Optional[bool]:
    """
        Instantiate curate objectt using command_dict from argument parser.

        :param command_dict:
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
                - flag: stores optional info for special data treatmen (i.e: ChEMBL input either a CHEMBL ID or a 
                list of CHEMBL ID in a file.)
    """
    
    # safety check if curation endpoint exists
    output_dir = utils.curation_tree_path(command_dict['endpoint'])
    if not os.path.isdir(output_dir):
        LOG.error("Endpoint name not found in model repository.\n")
        return
    
    # check of metadata

    if command_dict['metadata']:
        metadata_ = command_dict['metadata']
        if (command_dict['molecule_identifier'] in metadata_) or (command_dict['structure_column'] in metadata_):
            LOG.error("datacur curate : metadata can't contain the ID nor the SMILES column names.\n")
            return
    else:
        metadata_ = None
    
    # call of curation functions
    curation = metacur.Metacuration(data_input=command_dict['data_input'], 
                                    molecule_identifier=command_dict['molecule_identifier'],
                                    structure_column=command_dict['structure_column'],
                                    output_dir=output_dir,
                                    endpoint=command_dict['endpoint'],
                                    metadata=metadata_,
                                    separator=command_dict['separator'],
                                    remove_problematic=command_dict['remove_problematic'],
                                    curation_type=command_dict['curation_type'],
                                    flag=command_dict['flag'])
    
    if command_dict['curation_type'] == 'chem':
        curation.curate_data()
        curation.write_output_curation_data()
    elif command_dict['curation_type'] == 'htt':
        curation.chemical_curation()
        curation.write_chem_output_and_x_matrix()
        
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