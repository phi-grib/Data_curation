"""
    Different functions to handle miscelaneous issues

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 22/02/2021, 13:00 PM
"""

__modules__ = None

import os
import pathlib
import sys
import yaml

from typing import Union, Tuple

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