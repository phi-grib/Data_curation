"""
    Different functions to handle miscelaneous issues

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 22/02/2021, 13:00 PM
"""

__modules__ = None

import hashlib
import numpy as np
import os
import pathlib
import pickle
import random
import string
import sys
import time
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
        return False, 'unable to obtain configuration file'

    if conf['config_status']:
        items = ['curation_repository_path']
        for i in items:
            try:
                conf[i] = os.path.abspath(conf[i])
            except:
                return False, "Configuration file incorrect. Unable to convert {} to a valid path.".format(conf[i])

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