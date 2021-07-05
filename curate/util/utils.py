"""
    Different functions to handle miscelaneous issues

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 22/02/2021, 13:00 PM
"""

__modules__ = None

import os
import pandas as pd
import pathlib
import sys
import yaml

from rdkit import Chem
from rdkit.Chem import PandasTools
from typing import Union, Tuple, Optional

def isSingleThread() -> Union[str, bool]:
    """
        Checks if config file points only to a single thread.

        :return single_thread:
        :return False:
    """

    success, config = read_config()

    if success:
        if 'single_thread' in config:
            single_thread = config['single_thread']
            return single_thread
    
    return False

def config_test() -> bool:
    """ 
        Checks if datacur has been configured
        reading the config.yml and checking the config_status flag.
    """

    success, config = read_config()

    if success:
        if isinstance(config['config_status'], bool):
            if config['config_status']:
                if os.path.isdir(config['curation_repository_path']):
                    return
    
    sys.exit()

def read_config() -> Union[Tuple[bool, str], Tuple[bool, dict]]:
    """
        Reads configuration file "config.yaml" and checks
        sanity of model repository path.

        :return bool: returns True when config is read
        :return conf: config file
    """

    try:
        source_dir = os.path.dirname(os.path.dirname(__file__)) 
        config_nam = os.path.join(source_dir,'config.yaml')
        with open(config_nam,'r') as f:
            conf = yaml.safe_load(f)
    except Exception as e:
        return False, e

    if conf is None:
        return False, 'unable to obtain configuration file\n'

    if conf['config_status']:
        items = ['curation_repository_path']
        for i in items:
            try:
                conf[i] = os.path.abspath(conf[i])
            except:
                return False, "Configuration file incorrect. Unable to convert {} to a valid path.\n".format(conf[i])

    return True, conf

def write_config(config: dict):
    """
        Writes the configuration to disk.

        :param config: configuration passed as dictionary
    """

    config['config_status'] = True
    source_dir = os.path.dirname(os.path.dirname(__file__))

    with open(os.path.join(source_dir,'config.yaml'), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

def set_repositories(curation_path: str):
    """
        Set the curation repository path.
        This is the dir where datacur is going to create and load curations.

        :param curation_path:
    """
    
    source_dir = os.path.dirname(os.path.dirname(__file__))

    with open(os.path.join(source_dir,'config.yaml'), 'r') as f:
        configuration = yaml.safe_load(f)

    new_curation_path = pathlib.Path(curation_path)

    configuration['curation_repository_path'] = str(new_curation_path.resolve())
    
    write_config(configuration)

def curation_repository_path() -> str:
    """
        Returns the path to the root of the curation repository,
        containing all curations.

        :return cur_path: path to curation repository
    """

    success, config = read_config()
    cur_path = config['curation_repository_path']
    
    return cur_path

def curation_tree_path(curation_path: str) -> str:
    """
        Returns the path to the curation given as argument, containg all versions

        :param curation_path: directories to curation path endpoint

        :return path_to_return: complete tree path to specific curation endpoint
    """

    path_to_return = os.path.join(curation_repository_path(), curation_path)

    return path_to_return

def set_curation_repository(path: str = None):
    """
        Set the curation repository path.
        This is the dir where datacur is going to create and load curated datasets.
        if path is None, curation dir will be set to the default in the
        datacur root directory.

        :param path: string indicating the path of the curation repository
    """

    source_dir = os.path.dirname(os.path.dirname(__file__)) 
    with open(os.path.join(source_dir,'config.yaml'),'r') as f:
        configuration = yaml.safe_load(f)

    if path is None:  # set to default path
        model_root_path = os.path.join(
            pathlib.Path(__file__).resolve().parents[1],
            'curation/')
        configuration['curation_repository_path'] = str(model_root_path)
    else:
        new_path = pathlib.Path(path)
        configuration['curation_repository_path'] = str(new_path.resolve())

    write_config(configuration)

def format_output(data: pd.DataFrame, outfile_type: str, outfile_path: str, smiles_column: str = None, identifier: str = None):
        """
            Gives the desired format to the input data.
            Requires output file path and type of file (Excel, CSV, TSV, sdf, json)

            :param data: dataframe containing the data to be processed
            :param outfile_type: type of file to create
            :param outfile_path: output file path
            :param smiles_column: SMILES column in the dataframe to be processed (optional)
        """

        if 'sdf' in outfile_type.lower():
            write_sdf(data, outfile_path, smiles_column, identifier)
        elif 'xlsx' in outfile_type.lower() or 'excel' in outfile_type.lower():
            output_name_format = '.'.join([outfile_path,'xlsx'])
            data.to_excel(output_name_format)
        elif 'csv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_path,'csv'])
            data.to_csv(output_name_format, sep=',')
        elif 'tsv' in outfile_type.lower():
            output_name_format = '.'.join([outfile_path,'tsv'])
            data.to_csv(output_name_format, sep='\t')
        elif 'json' in outfile_type.lower():
            output_name_format = '.'.join([outfile_path,'json'])
            data.to_json(path_or_buf = output_name_format, orient = 'index')

def write_sdf(data: pd.DataFrame, outfile_name: str, smiles_column: str, identifier: str):
        """
            Prepares curated data to be converted into sdf file using
            PandasTools. Returns non processed molecules in excel format.

            :param data: Dataframe to be written
            :param smiles_column: SMILES column in the dataframe to be processed
            :param outfile_name: output file name
        """

        output_name_format = '.'.join([outfile_name,'sdf'])
        cur_data = prepare_data_for_sdf(data, smiles_column)
        
        PandasTools.WriteSDF(cur_data, output_name_format, molColName='ROMol', properties=list(cur_data.columns), idName=identifier)

def prepare_data_for_sdf(data: pd.DataFrame, smiles_column: str) -> Optional[pd.DataFrame]:
    """
        Prepares the data to be converted to sdf.

        :param data: Dataframe to be treated
        :param smiles_column: SMILES column in the dataframe to be processed

        :return data: dataframe with new columns added before being converted into sdf.
    """

    PandasTools.AddMoleculeColumnToFrame(data, smiles_column)
    no_mol = data[data['ROMol'].isna()]
    data.drop(no_mol.index, axis=0, inplace=True)
    data.loc[:,'ROMol'] = [Chem.AddHs(x) for x in data['ROMol'].values.tolist()]
    
    # if no_mol.empty is False:
    #     non_processed_path = '/'.join([output_dir,'Non_processed_molecules'])
    #     format_output(data = no_mol, outfile_type = 'xlsx', outfile_path = non_processed_path)

    return data